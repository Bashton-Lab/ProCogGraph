#!/usr/bin/env python

"""
Downloads and validates the external reference data files needed to run
the ProCogGraph Nextflow pipeline from scratch (see docs/installation.md
and nextflow/nextflow.config:85-103), driven by
reference_data_manifest.yaml. See docs/reference_data_download_plan.md
for the design rationale.

Example usage:
    python3 download_reference_data.py --data_dir /path/to/data_dir
    python3 download_reference_data.py --data_dir /path/to/data_dir --dry-run
    python3 download_reference_data.py --data_dir /path/to/data_dir --only enzyme_dat,rhea2ec
    python3 download_reference_data.py --data_dir /path/to/data_dir --force
"""

import argparse
import csv
import gzip
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
import yaml

DEFAULT_MANIFEST = Path(__file__).parent / "reference_data_manifest.yaml"
RETRY_ATTEMPTS = 3
RETRY_BACKOFF_SECONDS = 5
REQUEST_TIMEOUT_SECONDS = 120


def load_manifest(manifest_path):
    with open(manifest_path) as handle:
        manifest = yaml.safe_load(handle)
    return manifest["entries"]


def download_with_retry(url, dest_path):
    last_exception = None
    for attempt in range(1, RETRY_ATTEMPTS + 1):
        try:
            with requests.get(url, stream=True, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                response.raise_for_status()
                tmp_path = dest_path.with_suffix(dest_path.suffix + ".part")
                with open(tmp_path, "wb") as out_file:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        out_file.write(chunk)
                tmp_path.rename(dest_path)
            return
        except (requests.RequestException, OSError) as exc:
            last_exception = exc
            if attempt < RETRY_ATTEMPTS:
                print(f"    download failed (attempt {attempt}/{RETRY_ATTEMPTS}): {exc} - retrying in {RETRY_BACKOFF_SECONDS}s")
                time.sleep(RETRY_BACKOFF_SECONDS)
    raise RuntimeError(f"failed to download {url} after {RETRY_ATTEMPTS} attempts") from last_exception


def check_size(dest_path, min_size_bytes):
    if min_size_bytes is None:
        return
    actual_size = dest_path.stat().st_size
    if actual_size < min_size_bytes:
        raise RuntimeError(
            f"{dest_path.name} is {actual_size} bytes, below the expected minimum of "
            f"{min_size_bytes} bytes - likely a truncated download, an error page, or a "
            f"stale URL. Not treating this as a successful fetch."
        )


def post_process_gunzip(dest_path):
    gz_path = dest_path.with_name(dest_path.name + ".gz.tmp")
    dest_path.rename(gz_path)
    with gzip.open(gz_path, "rb") as src, open(dest_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    gz_path.unlink()


def post_process_csv_to_tsv_gz(dest_path):
    """SIFTS' only live source is comma-separated and uncompressed, but the
    pipeline (download_mmcif.py) reads it tab-separated from a .gz file."""
    csv_path = dest_path.with_name(dest_path.name + ".csv.tmp")
    dest_path.rename(csv_path)
    with open(csv_path, newline="") as src, gzip.open(dest_path, "wt", newline="") as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst, delimiter="\t")
        for row in reader:
            writer.writerow(row)
    csv_path.unlink()


def derive_chebi_results(data_dir, entry, manifest_by_name, force, skip_existing):
    """Reconstruct ChEBI_Results.tsv (ChEBI ID/name/SMILES/KEGG COMPOUND
    ACCESSION for ChEBI entries with a KEGG cross-reference and a known
    structure) from bulk ChEBI flat files, instead of ChEBI's advanced
    search UI. See docs/reference_data_download_plan.md #2.

    Schema notes (verified against live ChEBI flat_files, 2026-08 - these
    are NOT documented anywhere obvious, so don't "simplify" this without
    re-checking against the real files first):
      - all ChEBI flat file columns are lowercase.
      - database_accession.tsv.gz's `type` column is a mix of CAS/CITATION/
        MANUAL_X_REF/REGISTRY_NUMBER - it does NOT contain a KEGG-specific
        type value. KEGG COMPOUND cross-references are `type ==
        "MANUAL_X_REF"` rows whose `source_id` is 45 in source.tsv.gz (not
        downloaded separately here - confirmed by inspecting source.tsv.gz
        directly; a stable ChEBI-internal primary key, not expected to
        change). Only KEGG COMPOUND is used here, not KEGG GLYCAN (47) or
        KEGG DRUG (46), to match get_ec_information.py's actual merge key
        (nextflow/bin/get_ec_information.py:536, `right_on = "KEGG COMPOUND
        ACCESSION"` - the original manual export also had separate
        KEGG GLYCAN/KEGG DRUG columns, but neither is ever read).
      - structures.tsv.gz has a `smiles` column directly (one row per
        compound_id, not melted by structure type) - no TYPE=="SMILES"
        filtering needed or possible.
      - the "ID" column MUST be formatted as "CHEBI:<n>" (not a bare
        integer) - get_ec_information.py:543 sets `ligand_db` directly
        from this column, and :682 later does
        `ligand_db.str.extractall("CHEBI:([0-9]+)")` to feed ChEBI
        has_role-based cofactor tagging. A bare integer silently produces
        zero regex matches there rather than an error.
    """
    KEGG_COMPOUND_SOURCE_ID = 45

    dependency_names = ["chebi_names", "chebi_database_accession", "chebi_structures"]
    dependency_paths = {}
    for dep_name in dependency_names:
        dep_entry = manifest_by_name[dep_name]
        dep_path = data_dir / dep_entry["target_filename"]
        if not dep_path.exists():
            print(f"    dependency {dep_name} not present, fetching it first")
            fetch_entry(data_dir, dep_entry, manifest_by_name, force, skip_existing)
        dependency_paths[dep_name] = dep_path

    names_df = pd.read_csv(dependency_paths["chebi_names"], sep="\t", compression="gzip")
    names_df = names_df.groupby("compound_id").agg({"name": "first"}).reset_index()

    accession_df = pd.read_csv(dependency_paths["chebi_database_accession"], sep="\t", compression="gzip")
    kegg_accession_df = accession_df.loc[
        (accession_df["type"] == "MANUAL_X_REF") & (accession_df["source_id"] == KEGG_COMPOUND_SOURCE_ID)
    ][["compound_id", "accession_number"]].rename(columns={"accession_number": "KEGG COMPOUND ACCESSION"})

    structures_df = pd.read_csv(dependency_paths["chebi_structures"], sep="\t", compression="gzip")
    smiles_df = structures_df.dropna(subset=["smiles"])[["compound_id", "smiles"]].drop_duplicates(
        subset="compound_id", keep="first"
    ).rename(columns={"smiles": "SMILES"})

    merged = kegg_accession_df.merge(smiles_df, on="compound_id", how="inner")
    merged = merged.merge(names_df, on="compound_id", how="left")
    merged["ID"] = "CHEBI:" + merged["compound_id"].astype(str)
    merged = merged.rename(columns={"name": "NAME"})
    merged = merged[["ID", "NAME", "SMILES", "KEGG COMPOUND ACCESSION"]].drop_duplicates()

    dest_path = data_dir / entry["target_filename"]
    merged.to_csv(dest_path, sep="\t", index=False)
    return dest_path


DERIVED_FUNCTIONS = {
    "derive_chebi_results": derive_chebi_results,
}

POST_PROCESS_FUNCTIONS = {
    "gunzip": post_process_gunzip,
    "csv_to_tsv_gz": post_process_csv_to_tsv_gz,
}


def fetch_entry(data_dir, entry, manifest_by_name, force, skip_existing):
    name = entry["name"]
    source_type = entry["source_type"]
    dest_path = data_dir / entry["target_filename"] if entry.get("target_filename") else None

    if dest_path is not None and dest_path.exists() and skip_existing and not force:
        print(f"[skip] {name}: {dest_path.name} already exists")
        return "skipped", dest_path

    if source_type == "manual":
        print(f"[manual] {name}: no automated source available.")
        if entry.get("note"):
            print(f"    {entry['note'].strip()}")
        return "manual", None

    if source_type == "needs_code_update":
        print(f"[fetched, needs code update] {name}: downloading {dest_path.name}, but pipeline code isn't ready to consume it yet.")
        if entry.get("note"):
            print(f"    {entry['note'].strip()}")
        download_with_retry(entry["url"], dest_path)
        check_size(dest_path, entry.get("min_size_bytes"))
        return "needs_code_update", dest_path

    if source_type == "derived":
        print(f"[derive] {name}: deriving {dest_path.name}")
        derive_fn = DERIVED_FUNCTIONS[entry["post_process"]]
        derive_fn(data_dir, entry, manifest_by_name, force, skip_existing)
        check_size(dest_path, entry.get("min_size_bytes"))
        return "derived", dest_path

    if source_type == "direct_url":
        print(f"[fetch] {name}: downloading {entry['url']}")
        download_with_retry(entry["url"], dest_path)
        post_process = entry.get("post_process", "none")
        if post_process != "none":
            POST_PROCESS_FUNCTIONS[post_process](dest_path)
        check_size(dest_path, entry.get("min_size_bytes"))
        return "fetched", dest_path

    raise ValueError(f"unknown source_type '{source_type}' for entry '{name}'")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data_dir", required=True, type=str, help="output directory (matches params.data_dir)")
    parser.add_argument("--manifest", type=str, default=str(DEFAULT_MANIFEST), help="path to reference_data_manifest.yaml")
    parser.add_argument("--only", type=str, default=None, help="comma-separated entry names to fetch (default: all)")
    parser.add_argument("--force", action="store_true", help="re-download even if the target file already exists")
    parser.add_argument("--include-optional", action="store_true", help="also fetch entries marked optional (e.g. cofactor-coverage data, not required by the current pipeline)")
    parser.add_argument("--dry-run", action="store_true", help="print what would be fetched/skipped without downloading")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    entries = load_manifest(args.manifest)
    manifest_by_name = {entry["name"]: entry for entry in entries}

    if args.only:
        requested_names = set(args.only.split(","))
        unknown = requested_names - set(manifest_by_name)
        if unknown:
            print(f"error: unknown entry name(s): {', '.join(sorted(unknown))}", file=sys.stderr)
            sys.exit(1)
        entries = [entry for entry in entries if entry["name"] in requested_names]
    elif not args.include_optional:
        entries = [entry for entry in entries if not entry.get("optional", False)]

    if args.dry_run:
        print(f"Would process {len(entries)} entries into {data_dir}:\n")
        for entry in entries:
            target = entry.get("target_filename") or "(derived, no fixed target)"
            exists = (data_dir / entry["target_filename"]).exists() if entry.get("target_filename") else False
            action = "skip (exists)" if exists and not args.force else entry["source_type"]
            print(f"  [{action:20s}] {entry['name']:35s} -> {target}")
        return

    results = {}
    log_lines = [f"# download_reference_data.py run at {datetime.now(timezone.utc).isoformat()}"]
    for entry in entries:
        try:
            status, dest_path = fetch_entry(data_dir, entry, manifest_by_name, args.force, skip_existing=not args.force)
        except Exception as exc:
            print(f"[FAILED] {entry['name']}: {exc}")
            status, dest_path = "failed", None
        results.setdefault(status, []).append(entry["name"])
        source_note = entry.get("url") or "(derived/manual)"
        log_lines.append(f"{entry['name']}\t{status}\t{source_note}\t{dest_path or ''}")

    log_path = data_dir / "download_manifest.log"
    with open(log_path, "a") as log_file:
        log_file.write("\n".join(log_lines) + "\n")

    print("\nSummary:")
    for status in ["fetched", "derived", "skipped", "needs_code_update", "manual", "failed"]:
        names = results.get(status, [])
        if names:
            print(f"  {status}: {len(names)} ({', '.join(names)})")

    if results.get("failed"):
        sys.exit(1)


if __name__ == "__main__":
    main()
