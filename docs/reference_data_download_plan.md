# Plan: Automate External Reference Data Download

**Status (2026-08-02, branch `feature/reference-data-download`, off
`fix/pipeline-wiring-bugs`):** Implemented as
`nextflow/bin/reference_data_manifest.yaml` +
`nextflow/bin/download_reference_data.py`, with `direct_url`/`derived`/
`manual`/`needs_code_update` source types, tested end-to-end against live
sources (not just written and assumed correct) during implementation.
Notable deviations/fixes found only by actually running it against real
data, not by re-reading this plan:

- The `ChEBI_Results.tsv` derivation logic in this plan's §2 was wrong in
  three ways only caught by inspecting the real files: ChEBI flat-file
  columns are lowercase, not uppercase as assumed;
  `database_accession.tsv.gz`'s `type` column has no KEGG-specific value
  at all (KEGG cross-references are `type == "MANUAL_X_REF"` rows filtered
  by `source_id` against `source.tsv.gz`, a file this plan never
  mentioned); and `structures.tsv.gz` has a direct `smiles` column per
  compound, not a melted `TYPE`/`STRUCTURE` format. Also: the derived `ID`
  column must be `CHEBI:<n>`-prefixed and KEGG-COMPOUND-only (not also
  KEGG GLYCAN, which this plan's §2 originally included) to match what
  `get_ec_information.py` actually merges on and what its ChEBI
  `has_role`-based cofactor-tagging regex (`:682`) expects — a bare
  integer ID would have silently zeroed out that regex match rather than
  erroring. Verified by diffing derived output against the repo's existing
  manually-exported `data_files/ChEBI_Results.tsv`.
- SIFTS's only live source is a comma-separated, uncompressed `.csv`, not
  the tab-separated gzipped `.tsv.gz` `download_mmcif.py` reads — added a
  `csv_to_tsv_gz` post-process step (not anticipated by this plan) rather
  than changing the pipeline script.
- SCOP (`dir.cla.scop.1_75.txt`/`dir.des.scop.1_75.txt`), SCOP2
  (`scop2-cla-latest.txt`/`scop2-des-latest.txt`), the CCD
  (`components.cif.gz`), and the three RHEA TSVs could not be confirmed
  live during implementation (inconclusive search results / a directory
  listing rather than the file itself) — the manifest uses best-known
  conventional URLs for these, each marked `confidence: unverified` with a
  note. Run `--dry-run --only <name>` and then a real fetch for these
  specifically before depending on them; update the manifest if any 404s.
- `preprocess_rhea.py --rd_dir` needs a directory of Rhea `.rd` files with
  no bulk source identified — left as a `manual`-style gap
  (`rhea_rd_dir` entry), not solved here.
- The Pfam clan-file consolidation (`Pfam-A.clans.tsv.gz`) is downloaded
  (`needs_code_update` source type) but intentionally excluded from the
  script's "ready to run pipeline" framing, since `utils.get_pfam_annotations`
  doesn't parse the new format yet — that remains separate work, per this
  plan's original §2.

`docs/installation.md` updated to point at the script in place of the old
manual download table.

Currently, running the Nextflow pipeline from scratch requires manually
downloading ~20 external reference files (SIFTS, PDBe-KB, ExPASy, CATH,
SCOP, SCOP2, Pfam/InterPro, RHEA, ChEBI, CCD, PubChem) into a single
`data_dir`, per the table in [`docs/installation.md`](installation.md#L99).
This is tedious, error-prone (wrong filename/location silently breaks a
pipeline stage later), and undocumented in a machine-checkable way. This
plan proposes a single script that fetches and validates all of them.

## 1. Build a authoritative manifest of files needed

Cross-referencing the installation doc's table against actual
`nextflow.config` params (`nextflow/nextflow.config:86-102`) and script
argparse definitions turned up two gaps to fix as part of this work, not
after:

- **`pfamA.txt.gz`** is required by `params.pfam_a_file`
  (`nextflow.config:95`, consumed by `utils.get_pfam_annotations`), but is
  **not listed** in the installation doc's download table at all.
- **`scop2-des-latest.txt`** is listed in the doc's table but has **no
  corresponding `params.scop2_descriptions_file` entry in
  `nextflow.config`**, even though `process_all_pdb_contacts.py` and
  `produce_neo4j_files.py` both take `--scop2_descriptions_file` as a
  required argument (this is fix_plan.md item 2 — will need a
  `nextflow.config` param added regardless of this plan, but the download
  script must fetch the file either way).

So step 1 is: build one manifest (e.g.
`nextflow/bin/reference_data_manifest.yaml`) that is the single source of
truth for "what files does the pipeline need, and where do they come
from," covering all params referenced in `nextflow.config:86-102` plus
`scop2_descriptions_file`. The installation doc's table and
`nextflow.config` should both be generated from / checked against this
manifest, not maintained separately, so they can't drift again.

Each manifest entry needs: `name`, `target_filename`, `param_name` (the
matching `nextflow.config` param, for cross-checking), `source_type`
(`direct_url` / `latest_release_dir` / `manual`), `url`, and optionally
`post_process` (e.g. `gunzip`, `none`).

## 2. Handle three distinct source patterns

The current file list isn't uniformly "download this URL" — group by
pattern so the script can treat each correctly:

- **Direct, stable URL** (majority of files): SIFTS, PDBe-KB assemblies,
  RHEA files, ChEBI names/relations, CCD — just `curl`/`requests.get`.
- **"Latest release" directory listings**: CATH
  (`cath-classification-data/`), SCOP2 (`pdbe/scop/download`), Pfam/InterPro
  (`current_release/`) don't have one fixed filename/version — the script
  needs to resolve the actual file to fetch (e.g. parse a directory index
  or use a known-stable "latest" alias URL if the source provides one).
  Where the source has no stable "latest" alias, log the resolved
  version/date so re-runs are reproducible and it's visible when upstream
  data has moved on.

  **Path correction (checked live, 2026-08):** Pfam has also restructured.
  `clan_membership.txt.gz` and `clan.txt.gz` — both named explicitly in
  `docs/installation.md`'s download table — **no longer exist** at
  `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/`. Clan data is
  now consolidated into a single file, **`Pfam-A.clans.tsv.gz`** (541K).
  `pfamA.txt.gz` (13M, already flagged as missing from the doc's table in
  step 1) is still present under that name. The manifest entries for both
  `pfam_clan_rels`/`pfam_clans` params (`nextflow.config:93-94`) need to
  point at the one consolidated file instead of two, and
  `utils.get_pfam_annotations` (`nextflow/bin/utils.py:232-241`), which
  currently expects two separate files (`clan_membership_file`,
  `clan_info_file`), will need updating to parse the merged format — this
  is a pipeline-code change, not just a download-path fix, so it should be
  scoped alongside `fix_plan.md` rather than silently patched inside the
  download script.
- **Derived from bulk files** (not manual): `ChEBI_Results.tsv` was
  originally produced via ChEBI's advanced search UI
  (`docs/installation.md:119`), which looked non-scriptable at first
  glance. Tracing its actual consumer
  (`nextflow/bin/get_ec_information.py:531-536`) shows it only needs four
  columns — ChEBI ID, ChEBI name, SMILES, and KEGG COMPOUND
  ACCESSION — restricted to ChEBI entries that (a) have a KEGG
  COMPOUND/GLYCAN cross-reference and (b) have a structure on file. Both
  of those facts are published by ChEBI as plain bulk flat files.

  **Path correction (checked live, 2026-08):** ChEBI has since
  restructured its FTP layout entirely. The
  `Flat_file_tab_delimited/` directory referenced by
  `docs/installation.md` no longer exists — it's now
  `https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/`, and the
  filenames changed too:
  - `chebi_names.tsv.gz` → **`names.tsv.gz`** (9.0M compressed)
  - `relation.tsv` → **`relation.tsv.gz`** (now gzipped, 2.6M)
  - the cross-reference table needed here is **`database_accession.tsv.gz`**
    (now gzipped, 3.8M) — filter `TYPE` to `KEGG COMPOUND accession` /
    `KEGG GLYCAN accession`.
  - the structure table is **`structures.tsv.gz`** (not `.csv.gz` —
    renamed, and larger than assumed at 88M compressed) — filter `TYPE`
    to `SMILES`.

  Any manifest entry or doc text written before this check needs the
  corrected path/filename/extension, not just the derivation logic. This
  is exactly the kind of drift the manifest in step 1 is meant to catch
  early via a periodic `--dry-run` / size-sanity check, rather than
  silently downloading a 404 page or failing a `pd.read_csv`.

  So this isn't a `source_type: manual` entry — it's a `derived` entry:
  download both bulk files (now `.gz`), inner-join them on ChEBI compound
  ID, merge in the name from `names.tsv.gz`, and write out
  `ChEBI_Results.tsv` in the same shape the pipeline already expects. This
  is strictly better than the original manual export too: it's
  reproducible and pinned to a specific ChEBI release rather than a
  point-in-time UI query result.
- **Manual-only, no stable download URL**: none currently identified. Keep
  `source_type: manual` as a category in the manifest schema for future
  entries, but the manifest should end up with zero of them once the
  `derived` handling above is implemented. If a genuinely manual entry
  ever shows up, the script should skip it, print clear instructions, and
  continue rather than fail the whole run.

## 3. Script design

New script: `nextflow/bin/download_reference_data.py`

- `--data_dir` (required) — output directory, matches `params.data_dir`.
- `--manifest` (default: bundled `reference_data_manifest.yaml`).
- `--only NAME[,NAME...]` — fetch a subset (useful when only one upstream
  source has updated).
- `--skip-existing` (default true) / `--force` — re-download even if the
  target file already exists, for refreshing stale reference data.
- `--dry-run` — print what would be fetched/skipped without downloading.
- Behaviour per file:
  1. Skip if present and `--skip-existing` (default), unless `--force`.
  2. Download with retry/backoff (a handful of these hosts — EBI FTP,
     Zenodo-adjacent mirrors — are occasionally flaky) and a sane timeout.
  3. Verify the download is non-empty / matches expected minimum size as a
     sanity check (no fixed checksums available from most of these
     upstream sources, so this is a best-effort integrity check, not a
     hash).
  4. Apply `post_process` (e.g. leave `.gz` files as-is where the pipeline
     expects gzip, since several params like `sifts_file`/`interpro_xml`
     are consumed compressed). For `derived` entries, `post_process` is the
     join/filter logic described in step 2 (e.g. a small pandas function
     registered per entry) rather than a simple decompress.
  5. Log source URL + resolved version/date actually fetched, into a
     `download_manifest.log` written to `data_dir`, so a user can tell
     later exactly what version of each upstream file they're running
     against — useful for reproducing / debugging a specific pipeline run.
- At the end, print a summary table: fetched / derived / skipped /
  manual-required, and exit non-zero only if a required (non-manual) file
  failed.

## 3a. New source: PDBe-KB / PDBeChem v2 ligand-function data (for cofactor coverage)

This is *not* one of the ~20 files currently required by the pipeline —
it's a new source needed only once
[`docs/cofactor_coverage_plan.md`](cofactor_coverage_plan.md) is
implemented, added here so the manifest stays the single place tracking
every external file the pipeline (present or planned) depends on.

`docs/installation.md` already lists one PDBe-KB file
(`assemblies_data.csv`, from `pdbe-kb/complexes/`) but that's a different
resource — assembly preference data, unrelated to ligand function. The
relevant PDBe-KB/PDBeChem v2 resource for cofactor coverage is separate:
`https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/additional_data/pdb_ligand_interactions/`,
containing `interacting_chains_with_ligand_functions.tsv` (202M) and
`pdb_bound_molecules.tsv` (137M). These are bulk, structure/bound-entity
level outputs of PDBe's own **RelLig** pipeline, which classifies every
observed PDB ligand as reactant-like, cofactor-like, or drug-like using
the same PARITY method ProCogGraph already implements — see
`cofactor_coverage_plan.md` for why this is the primary source to use for
that work, ahead of the UniProt-REST approach originally proposed there.
- `source_type: direct_url`, two files, no per-accession chunking needed
  (unlike the UniProt REST approach this replaces as primary source).
- Given the file sizes (well over the small reference tables elsewhere in
  this manifest — 200M+ each), the size-sanity check described in the
  script design below (§3) should use a floor generous enough not to
  false-flag on minor version growth, but tight enough to catch a
  truncated/failed download.

## 4. Integration points

- `docs/installation.md` — replace the current manual download table
  (lines 99-124) with instructions to run the new script, keeping a
  smaller reference table of what it fetches (generated from the manifest,
  see step 1) for transparency. No manual-step instructions should be
  needed once `ChEBI_Results.tsv` is handled as a `derived` entry.
- Do **not** wire this into `main.nf` as a pipeline process — reference
  data changes far less often than pipeline runs, and folding it into the
  Nextflow DAG would force a re-check/re-download on every run. Keep it as
  a separate, explicitly-invoked setup step, analogous to
  `setup_docker_linux.sh` for the Docker path.
- Add the new script's dependencies (`requests`, `pyyaml` if not already
  present) to `nextflow/envs/procoggraph.yaml` if missing.

## 5. Verification

- `--dry-run` against the manifest should list exactly the files currently
  documented in `docs/installation.md`'s table, plus `pfamA.txt.gz`, minus
  `clan_membership.txt.gz`/`clan.txt.gz` (replaced by the single
  `Pfam-A.clans.tsv.gz`), and with `chebi_names.tsv.gz`/`relation.tsv`
  pointing at their corrected `flat_files/names.tsv.gz` /
  `flat_files/relation.tsv.gz` locations (see §2 above).
- A real run into a clean `data_dir` should produce a directory usable
  directly by `nextflow run main.nf` (module fix_plan items 1-2 aside),
  i.e. every `params.*_file`/`*_dir` referenced in
  `nextflow.config:86-102` (plus `scop2_descriptions_file` once added)
  resolves to an existing, non-empty file.
- Re-running with no flags should skip everything (fast, idempotent);
  re-running with `--force` should refresh everything and update
  `download_manifest.log`.

## Suggested order of work

1. Build the manifest (step 1) — this alone surfaces/documents the
   `pfamA.txt.gz` and `scop2_descriptions_file` gaps precisely, which is
   useful even before the script exists.
2. Implement the script for the direct-URL files first (majority of the
   list, lowest risk).
3. Add "latest release" resolution for CATH/SCOP2/Pfam/InterPro.
4. Add the `derived` handling for `ChEBI_Results.tsv` (fetch
   `database_accession.tsv` + `structures.csv.gz`, join, write output) and
   validate its output against a copy of the original manually-exported
   file if one is still available, to catch any subtle format/column
   differences before relying on it.
5. Update `docs/installation.md` to point at the script.

This plan only covers *downloading* the reference data — it doesn't touch
the pipeline wiring bugs already tracked in
[`docs/fix_plan.md`](fix_plan.md). Fixing fix_plan.md item 2
(`scop2_descriptions_file` missing from `nextflow.config`) can happen
independently, but should land before or alongside this work so the
manifest's `param_name` cross-check has something to check against.
