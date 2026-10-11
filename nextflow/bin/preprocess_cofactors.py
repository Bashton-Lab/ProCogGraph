#!/usr/bin/env python

"""
Preprocesses cofactor-EC association data into a cognate-ligand-shaped
dataframe, for the coverage-gap plan in docs/cofactor_coverage_plan.md.
Mirrors preprocess_rhea.py's pattern: a standalone script producing a
pickle (cofactor_ligands_df.pkl) that get_ec_information.py reads in and
concatenates alongside its existing reaction-derived sources, tagged
ligand_source="cofactor".

Two combined, complementary sources (see the plan doc for the coverage
analysis behind this choice):

  1. CoFactor DB 2010 + BRENDA, via RelLig's own vendored cofactor_ec.csv /
     cofactors_details.json (Apache-2.0, pdberellig project) - 27 organic
     cofactor classes, each EC exact (4-digit), resolved to a SMILES via
     its representative PDB CCD code.
  2. UniProt's `COFACTOR` annotation, bulk-pulled (structured
     Name=...;Xref=ChEBI:... rows only - free-text-only rows, e.g.
     chlorophyll's photosystem entries, are out of scope here and are a
     deferred separate LLM-extraction follow-on per the plan doc).

Both sources' EC values get classified by how many of their 4 segments are
fully resolved (see utils.classify_ec_completeness): exact (4) passes
through, subsubclass-level partials ("N.N.N.-", 3) get broadcast to every
matching terminal EC via utils.broadcast_subsubclass_ec, and anything
coarser (class/subclass-level, 1 or 2) is dropped and logged rather than
broadcast - broadcasting those was checked and found to produce nonsense
(e.g. a class-level wildcard matching literally every terminal EC in that
class). See docs/cofactor_coverage_plan.md ("The broadcast problem, and
its fix") for the full reasoning and the numbers behind this rule.

Example usage:
    python3 preprocess_cofactors.py \
        --cofactor_ec_csv cofactor_ec.csv \
        --cofactor_details_json cofactors_details.json \
        --ccd_cif ccd.cif \
        --uniprot_cofactor_tsv uniprot_cofactor_annotations.tsv \
        --chebi_structures chebi_structures.tsv.gz \
        --chebi_names chebi_names.tsv.gz \
        --ec_dat enzyme.dat \
        --enzyme_class_file enzclass.txt \
        --outdir /path/to/output/directory
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd
from gemmi import cif
from rdkit.Chem import PandasTools

from utils import (
    process_ec_records,
    classify_ec_completeness,
    broadcast_subsubclass_ec,
    get_chem_comp_descriptors,
)

CHEBI_XREF_PATTERN = re.compile(r"Xref=ChEBI:CHEBI:(\d+)")


def load_cofactor_db_table(cofactor_ec_csv_path, cofactor_details_json_path):
    """Returns a DataFrame with one row per (EC, cofactor class) pair from
    CoFactor DB 2010 + BRENDA, columns: ec, cofactor_id, representative_ccd,
    source. cofactor_ec.csv is EC_NO/COFACTOR_ID/SOURCE; cofactors_details.json
    is a list of {id, representative, template, threshold} - representative
    is the PDB CCD code used to resolve a SMILES for that cofactor class."""
    cofactor_ec = pd.read_csv(cofactor_ec_csv_path)
    cofactor_ec.columns = [c.strip().strip('"') for c in cofactor_ec.columns]
    cofactor_ec["EC_NO"] = cofactor_ec["EC_NO"].astype(str).str.strip().str.strip('"')
    cofactor_ec["COFACTOR_ID"] = cofactor_ec["COFACTOR_ID"].astype(int)

    with open(cofactor_details_json_path) as handle:
        details = json.load(handle)
    details_df = pd.DataFrame(details)[["id", "representative"]].rename(
        columns={"id": "COFACTOR_ID", "representative": "representative_ccd"}
    )

    merged = cofactor_ec.merge(details_df, on="COFACTOR_ID", how="inner")
    merged = merged.rename(columns={"EC_NO": "ec", "COFACTOR_ID": "cofactor_id", "SOURCE": "source"})
    return merged[["ec", "cofactor_id", "representative_ccd", "source"]].drop_duplicates()


def parse_uniprot_cofactor_df(raw_df):
    """Takes the bulk UniProt TSV (columns: Entry, EC number, Cofactor) as
    already-loaded a DataFrame, and returns one row per (ec, chebi_id) pair
    for rows with at least one structured `Xref=ChEBI:CHEBI:<n>` in the
    Cofactor field. Rows with only a free-text Note (no Xref) are dropped
    here - confirmed in the plan doc these are ~0.7% of the bulk pull
    (chlorophyll's photosystem annotations among them) and are out of scope
    for this mechanical path."""
    records = []
    for _, row in raw_df.iterrows():
        cofactor_text = row.get("Cofactor")
        if not isinstance(cofactor_text, str):
            continue
        chebi_ids = CHEBI_XREF_PATTERN.findall(cofactor_text)
        if not chebi_ids:
            continue
        ec_field = row.get("EC number")
        if not isinstance(ec_field, str) or not ec_field.strip():
            continue
        ecs = [e.strip() for e in ec_field.split(";") if e.strip()]
        for ec in ecs:
            for chebi_id in chebi_ids:
                records.append({"ec": ec, "chebi_id": int(chebi_id)})
    if not records:
        return pd.DataFrame(columns=["ec", "chebi_id"])
    return pd.DataFrame(records).drop_duplicates()


def classify_and_split_ec_rows(df, terminal_ec_list, ec_col="ec"):
    """Splits a dataframe by how resolved its EC values are
    (utils.classify_ec_completeness): exact (4) rows pass through
    unchanged; subsubclass-level (3) rows are exploded into one row per
    matching terminal EC (utils.broadcast_subsubclass_ec); class/subclass-
    level (1 or 2) rows are dropped, not broadcast. Returns (kept_df,
    dropped_df) - dropped_df is for logging/auditability, not silent
    discard (docs/cofactor_coverage_plan.md's stated requirement)."""
    terminal_ec_list = list(terminal_ec_list)
    kept_rows = []
    dropped_rows = []

    for _, row in df.iterrows():
        ec = row[ec_col]
        level = classify_ec_completeness(ec)
        if level == 4:
            kept_rows.append({**row.to_dict(), "entry": ec})
        elif level == 3:
            matches = broadcast_subsubclass_ec(ec, terminal_ec_list)
            for terminal_ec in matches:
                kept_rows.append({**row.to_dict(), "entry": terminal_ec})
            if not matches:
                dropped_rows.append({**row.to_dict(), "reason": "subsubclass-level, no matching terminal EC"})
        else:
            dropped_rows.append({**row.to_dict(), "reason": f"too coarse to broadcast safely (level {level})"})

    kept_df = pd.DataFrame(kept_rows).drop(columns=[ec_col], errors="ignore") if kept_rows else pd.DataFrame(columns=list(df.columns) + ["entry"]).drop(columns=[ec_col], errors="ignore")
    dropped_df = pd.DataFrame(dropped_rows) if dropped_rows else pd.DataFrame(columns=list(df.columns) + ["reason"])
    return kept_df, dropped_df


def resolve_chebi_smiles_map(chebi_ids, chebi_structures_df, chebi_names_df):
    """Resolves a list of bare integer ChEBI IDs to SMILES/name, from the
    full ChEBI bulk structures/names files (not the pipeline's existing
    ChEBI_Results.tsv, which is deliberately narrowed to KEGG-COMPOUND-
    cross-referenced entries only - most cofactor ChEBI IDs, especially
    metal ions, have no KEGG COMPOUND cross-reference at all and would be
    silently dropped if resolved that way)."""
    ids = pd.Series(sorted(set(chebi_ids)), name="compound_id")
    structures = chebi_structures_df.dropna(subset=["smiles"]).drop_duplicates(subset="compound_id", keep="first")
    names = chebi_names_df.groupby("compound_id").agg({"name": "first"}).reset_index()

    resolved = ids.to_frame().merge(structures[["compound_id", "smiles"]], on="compound_id", how="inner")
    resolved = resolved.merge(names, on="compound_id", how="left")
    return resolved


def build_cofactor_ligands_df(cofactor_db_df, uniprot_df, terminal_ec_list, ccd_doc, chebi_structures_df, chebi_names_df):
    """Orchestrates both sources into a single dataframe shaped like
    get_ec_information.py's other cognate-ligand source frames (entry,
    compound_id, compound_name, ROMol, ligand_db, compound_reaction,
    ligand_source), ready to concat directly into cognate_ligands_df.
    Returns (cofactor_ligands_df, dropped_ec_log_df)."""
    cofactor_db_kept, cofactor_db_dropped = classify_and_split_ec_rows(cofactor_db_df, terminal_ec_list)
    uniprot_kept, uniprot_dropped = classify_and_split_ec_rows(uniprot_df, terminal_ec_list)

    rows = []

    if not cofactor_db_kept.empty:
        ccd_codes = cofactor_db_kept["representative_ccd"].dropna().unique().tolist()
        ccd_smiles = get_chem_comp_descriptors(ccd_doc, ccd_codes)
        for _, row in cofactor_db_kept.iterrows():
            smiles = ccd_smiles.get(row["representative_ccd"])
            if smiles is None:
                continue
            rows.append({
                "entry": row["entry"],
                "compound_id": row["representative_ccd"],
                "compound_name": row["representative_ccd"],
                "smiles": smiles,
                "ligand_db": f"CofactorDB:{row['cofactor_id']}",
            })

    if not uniprot_kept.empty:
        chebi_map = resolve_chebi_smiles_map(uniprot_kept["chebi_id"].unique(), chebi_structures_df, chebi_names_df)
        chebi_map = chebi_map.set_index("compound_id")
        for _, row in uniprot_kept.iterrows():
            chebi_id = row["chebi_id"]
            if chebi_id not in chebi_map.index:
                continue
            resolved = chebi_map.loc[chebi_id]
            rows.append({
                "entry": row["entry"],
                "compound_id": f"CHEBI:{chebi_id}",
                "compound_name": resolved["name"] if pd.notna(resolved.get("name")) else f"CHEBI:{chebi_id}",
                "smiles": resolved["smiles"],
                "ligand_db": f"CHEBI:{chebi_id}",
            })

    dropped_df = pd.concat([cofactor_db_dropped, uniprot_dropped], ignore_index=True)

    if not rows:
        empty = pd.DataFrame(columns=["entry", "compound_id", "compound_name", "ROMol", "ligand_db", "compound_reaction", "ligand_source"])
        return empty, dropped_df

    cofactor_ligands_df = pd.DataFrame(rows).drop_duplicates(subset=["entry", "compound_id", "smiles"])
    PandasTools.AddMoleculeColumnToFrame(cofactor_ligands_df, smilesCol="smiles", molCol="ROMol")
    cofactor_ligands_df = cofactor_ligands_df.loc[cofactor_ligands_df["ROMol"].notna()].copy()
    cofactor_ligands_df["compound_reaction"] = ""
    cofactor_ligands_df["ligand_source"] = "cofactor"
    cofactor_ligands_df = cofactor_ligands_df[
        ["entry", "compound_id", "compound_name", "ROMol", "ligand_db", "compound_reaction", "ligand_source"]
    ].reset_index(drop=True)

    return cofactor_ligands_df, dropped_df


def main():
    parser = argparse.ArgumentParser(description="Preprocess cofactor-EC association data")
    parser.add_argument("--cofactor_ec_csv", required=True, help="RelLig's vendored cofactor_ec.csv (EC_NO/COFACTOR_ID/SOURCE)")
    parser.add_argument("--cofactor_details_json", required=True, help="RelLig's vendored cofactors_details.json (cofactor class -> representative CCD code)")
    parser.add_argument("--ccd_cif", required=True, help="cif file containing the chemical component dictionary in mmcif format")
    parser.add_argument("--uniprot_cofactor_tsv", required=True, help="Bulk UniProt TSV (Entry, EC number, Cofactor columns)")
    parser.add_argument("--chebi_structures", required=True, help="ChEBI bulk structures.tsv.gz file")
    parser.add_argument("--chebi_names", required=True, help="ChEBI bulk names.tsv.gz file")
    parser.add_argument("--ec_dat", required=True, help="Path to enzyme.dat file from EXPASY")
    parser.add_argument("--enzyme_class_file", required=True, help="Path to enzyme_class file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    ec_records_df = process_ec_records(args.ec_dat, args.enzyme_class_file)
    terminal_ec_list = ec_records_df.TRANSFER.unique().tolist()

    cofactor_db_df = load_cofactor_db_table(args.cofactor_ec_csv, args.cofactor_details_json)

    uniprot_raw_df = pd.read_csv(args.uniprot_cofactor_tsv, sep="\t")
    uniprot_df = parse_uniprot_cofactor_df(uniprot_raw_df)

    ccd_doc = cif.read(args.ccd_cif)
    chebi_structures_df = pd.read_csv(args.chebi_structures, sep="\t", compression="gzip")
    chebi_names_df = pd.read_csv(args.chebi_names, sep="\t", compression="gzip")

    cofactor_ligands_df, dropped_ec_log = build_cofactor_ligands_df(
        cofactor_db_df, uniprot_df, terminal_ec_list, ccd_doc, chebi_structures_df, chebi_names_df
    )

    cofactor_ligands_df.to_pickle(f"{args.outdir}/cofactor_ligands_df.pkl")
    dropped_ec_log.to_csv(f"{args.outdir}/cofactor_dropped_ec_log.tsv", sep="\t", index=False)

    print(f"Cofactor ligands: {len(cofactor_ligands_df)} rows, {cofactor_ligands_df['entry'].nunique()} distinct ECs")
    print(f"Dropped EC rows (too coarse to broadcast, or unresolved): {len(dropped_ec_log)} - see cofactor_dropped_ec_log.tsv")


if __name__ == "__main__":
    main()
