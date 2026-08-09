# Fix Plan: Pipeline & Import Errors

Tracking doc for the errors identified in the first-pass review (2026-08-01).
Each item lists the defect, the fix, and how to verify it. Work top to
bottom — later Nextflow stages depend on earlier ones, so the pipeline
can't be smoke-tested end-to-end until items 1-3 are resolved.

**Status (2026-08-02, branch `fix/pipeline-wiring-bugs`):** Items 1-7
implemented. Item 8 (dead code) left untouched — needs an author decision,
not a mechanical fix.

**Additional bugs found during implementation, beyond this plan's original
scope**, all in `PROCESS_ALL_CONTACTS` and `PRODUCE_NEO4J_FILES` — the
same two processes items 1-3 were already fixing. Verified by diffing each
process's declared Nextflow `output:`/script CLI flags against what the
underlying Python script's argparse and `to_csv()` calls actually use:

- `PROCESS_ALL_CONTACTS` script call used `--contacts` but
  `process_all_pdb_contacts.py`'s argparse expects `--contacts_file` —
  would have failed immediately with an "unrecognized arguments" error.
- `PROCESS_ALL_CONTACTS` never passed `--scop2_descriptions_file` at all
  (same missing-param class of bug as item 2, confirmed to also apply
  here, not just to the file being absent from `nextflow.config`).
- `PROCESS_ALL_CONTACTS`'s declared outputs `scop2b_fa_pdb_residue_interactions.csv.gz`
  / `scop2b_sf_pdb_residue_interactions.csv.gz` don't match what the script
  actually writes — `scop2_fa_...` / `scop2_sf_...` (no "b") — a second,
  previously-uncaught instance of the same bug class as item 3.
- `PRODUCE_NEO4J_FILES` input `superfamily_domain_ownership` was declared
  as a process input but never actually passed as a `--superfamily_domain_ownership`
  flag on the script's command line at all.
- `PRODUCE_NEO4J_FILES`'s entire declared `output:` block (58 entries) was
  stale relative to the script: every filename used `.csv.gz` while the
  script writes `.tsv.gz` throughout, several files were renamed
  (`g3dsa_domains_nodes` → `gene3d_domains_nodes`, `scop2b_fa`/`scop2b_sf`
  → `scop2_fa`/`scop2_sf`, etc.), and 13 files the script writes
  (SCOP2 hierarchy nodes/rels, `superfamily_fold_rels`,
  `cath_topology_domain_rels`) weren't declared as outputs at all. The
  corrected 71-entry list was verified two ways: (1) extracting every
  `to_csv(f"...")` call from `produce_neo4j_files.py` directly, and (2)
  cross-checking against `nextflow/bin/import_neo4j_data.sh`, which already
  references the correct actual filenames (e.g. `gene3d_domains_nodes.tsv.gz`),
  confirming the import script was kept in sync with the Python script
  while `main.nf`'s output declarations were not.
- `PRODUCE_NEO4J_FILES` was also missing `scop2_domains_info`/
  `scop2_domains_description` inputs entirely (same root cause as item 2 —
  `produce_neo4j_files.py` requires `--scop2_domains_info_file` and
  `--scop2_descriptions_file`, neither of which existed as inputs or was
  passed on the command line).

None of this changes the plan's original fixes — it's the same items 1-3,
just more extensively broken than the original review (a static read,
without cross-referencing every script's argparse/output against
`main.nf`) caught. Fixed on the same branch as part of the same items.

## 1. Undefined `interpro_domain_ownership` in `PRODUCE_NEO4J_FILES`

**Status: fixed on `fix/pipeline-wiring-bugs`.**

- **Where:** [`nextflow/main.nf:205`](../nextflow/main.nf#L205) (call site),
  input block [`nextflow/main.nf:126-142`](../nextflow/main.nf#L126)
- **Problem:** the process call references
  `--interpro_domain_ownership ${interpro_domain_ownership}`, but no such
  variable is declared as a process input. Declared inputs are
  `superfamily_domain_ownership`, `g3dsa_domain_ownership`,
  `scop2b_sf_domain_ownership`, `scop2b_fa_domain_ownership`.
- **Fix:**
  1. Check `produce_neo4j_files.py`'s argparse definition to confirm the
     full, current set of expected `--*_domain_ownership` flags.
  2. Add any missing `path` inputs to the `PRODUCE_NEO4J_FILES` input block
     and wire them from the correct upstream channel in the workflow block
     (likely the CATH-Gene3D/SCOP2/InterPro domain ownership channels
     produced earlier in `main.nf`).
  3. Fix the call line to reference only variables that exist as inputs,
     with names matching what the script expects
     (`--superfamily_domain_ownership`, `--gene3dsa_domain_ownership`,
     `--scop2_sf_domain_ownership`, `--scop2_fa_domain_ownership`, etc.).
- **Verify:** `nextflow run main.nf -profile standard -stub-run` (or a real
  run against a small structure set) reaches `PRODUCE_NEO4J_FILES` without
  an undefined-variable/null error, and all expected `*_domain_ownership`
  files are non-empty in the output.

## 2. Missing SCOP2 descriptions file for `PROCESS_ALL_CONTACTS`

**Status: fixed on `fix/pipeline-wiring-bugs`.**

- **Where:** `PROCESS_ALL_CONTACTS` process call in `nextflow/main.nf`
  (around line 214); consumed by
  [`nextflow/bin/process_all_pdb_contacts.py`](../nextflow/bin/process_all_pdb_contacts.py)
  via `utils.get_scop2_domains_info`.
- **Problem:** the script requires a SCOP2 descriptions file path, but
  `main.nf` never passes one, and no corresponding `params.*` entry exists
  in `nextflow.config`.
- **Fix:**
  1. Determine where the SCOP2 descriptions file is expected to come from
     (a downloaded reference file, similar to other SCOP2/CATH reference
     data already staged elsewhere in the pipeline).
  2. Add a `params.scop2_descriptions_file` (or equivalent) entry to
     `nextflow.config`, wire it as a process input, and pass it on the
     `process_all_pdb_contacts.py` command line.
  3. Confirm the flag name matches the script's argparse definition exactly.
- **Verify:** `PROCESS_ALL_CONTACTS` runs without the
  `TypeError: expected str, bytes or os.PathLike object, not NoneType`
  currently produced, and SCOP2 domain info appears correctly in its output.

## 3. Output filename mismatch: `g3dsa_pdb_residue_interactions.csv.gz`

**Status: fixed on `fix/pipeline-wiring-bugs`.**

- **Where:** declared output
  [`nextflow/main.nf:91`](../nextflow/main.nf#L91) (`g3dsa_...`) vs. actual
  file written in
  [`nextflow/bin/process_all_pdb_contacts.py:343,346`](../nextflow/bin/process_all_pdb_contacts.py#L343)
  (`gene3dsa_...`).
- **Problem:** Nextflow will report a missing output file because the
  process declares a different filename than the script produces.
- **Fix:** Pick one name and use it consistently — recommend standardizing
  on `gene3dsa_...` (matches the script and is the less ambiguous
  abbreviation) and updating `main.nf:91` to match, rather than changing
  the Python script. Grep the rest of `main.nf` and `nextflow/bin/` for
  `g3dsa_pdb_residue_interactions` to make sure no other reference still
  uses the old name.
- **Verify:** the process's declared `output:` block and the file(s)
  actually written by the script are identical strings; pipeline proceeds
  past this process without a missing-output error.

## 4. Missing comma silently drops SCOP2/Pfam domains from filter

**Status: fixed on `fix/pipeline-wiring-bugs`.**

- **Where:** [`nextflow/bin/process_pdb_contacts.py:188`](../nextflow/bin/process_pdb_contacts.py#L188)
  ```python
  .isin(["CATH", "SCOP", "SCOP2B", "SCOP2" "Pfam", "InterPro"])
  ```
- **Problem:** adjacent string literals `"SCOP2" "Pfam"` concatenate into a
  single string `"SCOP2Pfam"` that matches nothing in the data. Structures
  with genuine SCOP2 or Pfam domain annotations are silently excluded from
  this filter — no error, just missing data downstream.
- **Fix:** add the missing comma:
  ```python
  .isin(["CATH", "SCOP", "SCOP2B", "SCOP2", "Pfam", "InterPro"])
  ```
- **Verify:** re-run on a structure known to have SCOP2 and/or Pfam
  annotations (spot-check via SIFTS) and confirm `mmcif_domains` /
  `domain_info_dataframe_filtered` now includes them. Since this bug
  produces no exception, add a quick assertion or log line counting rows
  per `xref_db` value so future regressions of this kind are caught
  immediately rather than silently.

## 5. Copy-paste error: SCOP class node mislabeled in Neo4j import

**Status: fixed on `fix/pipeline-wiring-bugs`.**

- **Where:** [`nextflow/bin/import_neo4j_data.sh:35`](../nextflow/bin/import_neo4j_data.sh#L35)
  ```
  --relationships=IS_IN_SCOP_CLASS=import/scop_fold_class_rels.tsv.gz \
  --nodes=IS_IN_SCOP_CLASS=import/scop_class_nodes.tsv.gz \
  ```
- **Problem:** the `--nodes` flag reuses the relationship type
  `IS_IN_SCOP_CLASS` as a node label instead of a proper label like
  `scopClass`. Every SCOP class node imported this way is mislabeled,
  breaking any Cypher query, constraint, or dashboard panel that expects
  the correct label.
- **Fix:** change the node line's label to `scopClass` (matching the
  naming convention of the sibling lines, e.g. `scopSuperfamily`,
  `scopFold`, `cathClass`). Double-check `import_neo4j_data_minimal.sh` and
  any Cypher/NeoDash queries under `procogdash/` don't already assume the
  broken label (i.e. that nothing was built to work around this bug).
- **Verify:** after a fresh `compose_build.yaml` import, run
  `MATCH (n:scopClass) RETURN count(n)` in Neo4j Browser and confirm it
  returns the expected non-zero count, and `MATCH (n:IS_IN_SCOP_CLASS)`
  returns zero.

## 6. Stale Zenodo version reference in docs (lower confidence)

**Status: fixed on `fix/pipeline-wiring-bugs`.**

- **Where:** `docs/installation.md` vs. the Zenodo record/version
  downloaded by `setup_docker_linux.sh` / `setup_docker_windows.ps1`.
- **Problem:** the doc appears to reference an older data release than what
  the setup scripts actually pull.
- **Fix:** confirm the current Zenodo record ID/version in the setup
  scripts, and update `docs/installation.md` (and any other doc/README
  mentioning a version) to match. Consider referencing the version from a
  single place (e.g. an env var or top-of-script constant) to avoid drift
  in future.
- **Verify:** every doc/script mentioning the Zenodo version string agrees.

## 7. Hardcoded author-specific conda path (lower confidence)

**Status: fixed on `fix/pipeline-wiring-bugs`.**

- **Where:** [`nextflow/main.nf:25`](../nextflow/main.nf#L25) vs. the
  portable `$baseDir`-relative conda path pattern used in
  `nextflow.config`.
- **Problem:** an absolute, machine-specific conda environment path breaks
  the pipeline for anyone who isn't the original author running it
  directly (outside the Docker-based prebuilt-flat-files path).
- **Fix:** replace the hardcoded path with the same `$baseDir`/relative
  pattern used elsewhere, so a fresh clone + `envs/procoggraph.yaml` setup
  works without editing `main.nf`.
- **Verify:** clone into a clean directory (or `git worktree add`), build
  the conda env from `envs/procoggraph.yaml`, and confirm `main.nf` doesn't
  reference any path outside the checkout.

## 8. Possible dead code (lower confidence, needs a decision)

- **Where:** `nextflow/bin/produce_minimal_neo4j_files.py`,
  `nextflow/bin/produce_flat_files.py`,
  `nextflow/bin/import_neo4j_data_minimal.sh`.
- **Problem:** none of these appear referenced from `main.nf`, the setup
  scripts, or docs. Either they're dead code left over from a prior
  approach, or they support an undocumented "minimal" install path.
- **Fix (needs author decision, not a mechanical fix):** determine which
  case it is.
  - If dead: remove the files (git history preserves them if needed later).
  - If it's a real, supported minimal-install path: wire it into `main.nf`
    or the setup scripts explicitly, and document it in
    `docs/installation.md`.
- **Verify:** N/A until the decision is made; then verify per whichever
  branch is chosen above.

## Suggested order of work

1. Items 1-3 (pipeline wiring) — the pipeline cannot run end-to-end without
   these, so fix and smoke-test together against a small structure set.
2. Item 4 (silent data-correctness bug) — high impact, cheap fix.
3. Item 5 (Neo4j import mislabeling) — verify against a fresh
   `compose_build.yaml` import once 1-3 are stable, since a full pipeline
   run + import is the only real way to confirm it.
4. Items 6-8 — low risk, can be done independently at any point; item 8
   needs a decision from the author before any file is touched.

No automated test suite exists yet. Once items 1-4 are fixed, consider
adding at least one smoke test (a tiny fixture structure run through
`process_pdb_contacts.py` / `process_all_pdb_contacts.py`) so regressions
like #4 fail loudly instead of silently next time.
