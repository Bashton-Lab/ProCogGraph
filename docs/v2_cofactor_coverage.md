# Cofactor Coverage in ProCogGraph v2

## Summary

ProCogGraph v1 derived cognate ligands exclusively from reaction
equations: an EC number's KEGG/Rhea reaction record supplies its
reactant/product compounds, which become candidate cognate ligands for
that EC. This is structurally blind to cofactors that participate
**catalytically rather than stoichiometrically** — prosthetic groups,
electron/light carriers, and structural metal centers that are never
written as a reactant or product of an EC's net reaction, and so can
never be discovered via that route no matter how complete the underlying
reaction database is.

ProCogGraph v2 adds a second, independent cognate-ligand source: a
combined EC→cofactor mapping built from the CoFactor database (Fischer,
Holliday & Thornton, *Bioinformatics* 2010) plus BRENDA, and from
UniProt's structured `COFACTOR` annotation. This is a strict addition
alongside the existing reaction-derived path, not a modification of it —
every cognate ligand row is now tagged with its provenance (`reaction`,
`cofactor`, or both, where a compound was independently discovered by
each route).

Benchmarked against the last pre-v2 build, this increases EC-ligand
coverage from 41,785 to 78,163 pairs (+87%) and the number of ECs with at
least one cognate ligand from 6,214 to 7,316 (+1,102 net), while leaving
the pre-existing reaction-derived coverage almost entirely intact (see
Results).

## Motivation

This is fine for substrates, products, and cofactors that *are*
stoichiometric participants (e.g. `NAD+ + substrate ⇌ NADH + product` —
NAD is written into the equation, so it's captured by v1 already). It
fails for the catalytic case: a photosystem's chlorophyll, for instance,
is never consumed or produced by the photosystem protein's own EC
reaction, so it is invisible to a reaction-equation-only pipeline
regardless of how good the underlying reaction data is.

## Method

### Sources

Two combined, independently-sourced EC→cofactor mappings:

1. **CoFactor database (2010) + BRENDA**, via the small vendored
   `cofactor_ec.csv`/`cofactors_details.json` tables from PDBe's RelLig
   project (Apache-2.0) — 27 organic cofactor classes (NAD, FAD, PLP, CoA,
   biotin, B12, heme A, SAM, molybdopterin, and others), each mapped to a
   representative PDB chemical-component code and, from there, a SMILES.
2. **UniProt's `COFACTOR` annotation**, bulk-pulled (not per-accession —
   a single `/uniprotkb/stream` query against all reviewed entries with
   both an EC number and a cofactor annotation), restricted to rows with
   a structured `Xref=ChEBI:CHEBI:<n>` (≈99.3% of the pull). This source
   is broader in chemical identity (109 distinct ChEBI cofactor
   identities vs. CoFactor DB's 27 classes — notably including metal
   centers such as `[4Fe-4S]` clusters) and catches EC/cofactor
   associations CoFactor DB's fixed, 2010-dated scope does not.

Measured real overlap between the two sources: only 1,201 of a combined
4,843 unique ECs are covered by both — genuinely complementary, not
redundant (concrete example: EC 1.1.1.10, a textbook NADP-dependent
enzyme, is covered by CoFactor DB but has zero UniProt cofactor
annotation on any reviewed entry).

### EC completeness and the broadcast rule

Not every source EC value is fully resolved to 4 digits. Each is
classified by how many segments are resolved before the first wildcard:

- **Exact (4/4)** — used directly.
- **Subsubclass-level (`N.N.N.-`)** — broadcast to every real terminal EC
  sharing that prefix (enzymes sharing all three leading digits generally
  do share cofactor chemistry).
- **Subclass- or class-level (`N.N.-.-`, `N.-.-.-`)** — **dropped, not
  broadcast**. An earlier attempt at broadcasting these produced a
  clearly wrong result (100% terminal-EC "coverage") by matching a
  class-level wildcard against every terminal EC in that class — e.g.
  claiming catalase and an unrelated NAD-dependent dehydrogenase share a
  cofactor purely because both are oxidoreductases. There is no
  chemically defensible way to narrow a class-level wildcard down to
  specific terminal ECs, so these are excluded rather than guessed at.

With this rule: 62.0% terminal-EC coverage from exact matches alone,
96.2% including safe subsubclass-level broadcast.

### Implementation

`nextflow/bin/preprocess_cofactors.py` builds the combined table
(mirroring `preprocess_rhea.py`'s standalone-script pattern) into
`cofactor_ligands_df.pkl`, which `get_ec_information.py
--cofactor_ligands` concatenates into `cognate_ligands_df` alongside the
existing Rhea/KEGG/ChEBI/PubChem/GlyTouCan sources, adding a
`ligand_source` column (`reaction` / `cofactor`, unioned where both
routes independently find the same compound). No existing matching code
(`get_pdb_parity.py`'s EC-keyed join) required modification — cofactor
rows are ordinary `cognate_ligands_df` rows, keyed by real terminal EC,
indistinguishable in shape from reaction-derived rows except for the
provenance tag.

## Results

Full end-to-end build against live current data (2026-08), compared
against the last pre-v2 `cognate_ligands_df.pkl` (dated 2024-07-18;
checksum-identical to the file that produced the published
[Zenodo v1-0-2 flat files](https://zenodo.org/records/14046116) and the
copy vendored into [AlphaCognate](https://github.com/m-crown/AlphaCognate)'s
`data/procoggraph_data/`):

| Metric | v1 (2024-07) | v2 (2026-08) | Change |
|---|---|---|---|
| EC–ligand pairs | 41,785 | 78,163 | **+87%** |
| ECs with ≥1 cognate ligand | 6,214 | 7,316 | **+1,102 net** (+1,111 gained / −9 lost) |
| Distinct ligand structures | 8,589 | 8,811 | +222 |

Of the 78,163 v2 pairs: 40,282 were already present in v1 (unchanged),
37,881 are new — of which 33,398 are corroborated by *both* the reaction
and cofactor paths independently finding the same compound, 2,406 are
reaction-path-only gains (unrelated to this work — normal upstream
Rhea/KEGG growth over the ~2-year gap between builds), and **2,077 exist
only because of the new cofactor path** — coverage that did not exist in
ProCogGraph v1 at all.

As a coherence check: the pre-existing, independent ChEBI `has_role`
cofactor-labelling mechanism (`isCofactor` column) — unrelated to this
work, present since v1 — shows a 6.3x increase in rows labelled
`Cofactor` (5,606 → 35,122). This is expected rather than circular: the
new cofactor-sourced rows are, by construction, literal cofactor
molecules (NAD, FAD, heme, etc.), and those largely already carry
ChEBI's own independent `has_role: cofactor` annotation — two separately-
sourced signals reinforcing each other.

### What didn't change (verification)

- Every pre-existing v1 EC-ligand pair's `ligand_source` backfills to
  `"reaction"` — the addition is strict, not a rewrite.
- `cognate_ligands_df.entry` contains zero wildcard/partial EC values
  after the broadcast step (hard invariant the downstream
  `get_pdb_parity.py` EC join depends on).
- Spot-checked cases match known biology: EC 1.1.1.10 gains
  FAD/TPP/NAD/NADP+/Mg²⁺/Zn²⁺; EC 1.97.1.12 (Photosystem I) gains a
  `[4Fe-4S] cluster` via the UniProt path; heme/heme b appear across the
  expected heme-dependent EC set.

### Known gaps (explained, not hidden)

**1,503 EC-ligand pairs (3.6% of the v1 total) present in v1 are absent
from v2.** Categorised (554 distinct lost structures):

- 180 involve wildcard-substituent SMILES (`*`) — partial/generic
  structures such as `[acyl-carrier protein]`-linked intermediates.
- 104 involve charge-separated porphyrin/macrocycle SMILES (`[N+]`/`[Mg`)
  — chlorophyll/bacteriochlorophyll/heme-biosynthesis-pathway
  intermediates specifically. **Root-caused 2026-10-11**: NOT an RDKit
  sanitization failure — every one of these 104 SMILES strings parses
  and sanitizes cleanly under current RDKit, ruling that out directly.
  Of the 104: 55 (53%) are a false positive of the pair-diff script
  itself — the identical compound (same name, same EC, same underlying
  Rhea reaction IDs, e.g. `RHEA:48765`/`RHEA:19734` for EC 1.1.1.396)
  is still present in v2, just written out with a different SMILES
  string for the same charge/coordination state (Rhea itself gives the
  same chlorophyll-pathway molecule different protonation/Mg-chelation
  notation across its own reaction entries — a pre-existing property of
  the upstream data, not something this work introduced). The script
  diffs on exact SMILES string equality, so these register as "lost"
  pairs when they're really the same ligand under a different textual
  representation. The remaining 49 (47%) are genuine compound-level
  gaps, but in every case the EC itself is still present in v2 with
  other ligands (cofactor and/or other reaction-derived) — no EC is
  lost because of this category. Consistent with the same normal
  upstream Rhea-revision pattern as the "other" category below, not a
  defect in this work.
- 270 are ordinary, otherwise-well-resolved compounds missing only for
  specific ECs — consistent with normal upstream Rhea/KEGG reaction-
  equation revision over the ~2-year gap (confirmed example: NADH is
  missing specifically for EC 1.1.1.96, but still resolves correctly for
  630 other ECs — the reaction equation for that one EC now cites NAD
  where it previously cited NADH).
- 9 ECs lost entirely (`1.14.14.140, 2.1.1.86, 2.4.1.129, 2.7.11.27,
  3.1.8.2, 4.2.1.78, 4.3.3.2, 4.3.3.3, 4.3.3.4`). **Root-caused
  2026-10-11**: not a real loss — all 9 are EC Nomenclature Committee
  "Transferred entry" numbers in the current `enzyme.dat` (confirmed:
  each one's `DE` line reads `Transferred entry: <new EC>`, e.g.
  `1.14.14.140` → `1.14.14.162`, `4.3.3.2` → `3.5.99.13`). The old
  baseline (2024-07) still carried ligand data directly under the old,
  now-retired EC numbers; v2 correctly follows the current EC
  nomenclature instead. Checked all 9 transfer targets directly: the
  exact same compounds (by name — e.g. dopamine, tryptamine, secologanin/
  strictosidine, water) are present under the new EC number in v2,
  usually with *more* ligands attached (the new EC also picks up
  cofactor-path compounds like heme, FAD, FMN, `[4Fe-4S]`). The pair-diff
  script missed this because `uniqueID` is a per-build index, not a
  stable cross-build compound identifier, so a naive join can't see the
  migration. Not a bug and nothing to fix — just a reminder that
  `enzyme.dat` transfers need resolving before any future old-vs-new
  comparison at the EC level.

**Chlorophyll as a photosystem cofactor is still not captured**, and this
is expected, not a bug: neither combined source resolves it mechanically.
CoFactor DB's 27 classes don't include chlorophyll at all, and UniProt's
photosystem cofactor annotations (e.g. `P56766`, EC 1.97.1.12) are
free-text `Note=` only, with no structured `Xref=ChEBI:...` — confirmed
directly against a live UniProt query (`cc_cofactor:"ChEBI:CHEBI:18230"`
returns zero reviewed entries). Chlorophyll's binding is described
qualitatively (variable stoichiometry, mixed a/a′ forms) rather than as a
clean 1:1 identity in both sources checked, which is plausibly *why*
neither gives it a structured entry. This affects a small, identified
long tail (742 of 111,819 UniProt cofactor+EC rows, ≈0.7%, are free-text
only) and is deferred to a separate LLM-assisted extraction pass, not
attempted mechanically here.

## Small-scale PARITY validation (2026-08-30/31)

The comparison in Results is at the `cognate_ligands_df.pkl` candidate-set
level (old vs. new cognate ligand lists) — it does not by itself confirm
that a real PDB structure's bound cofactor actually scores a PARITY match.
A small, hand-picked, three-structure run of the full downstream pipeline
(`download_mmcif.py` → `process_pdb_structure.py` → PDBe-Arpeggio →
`process_pdb_contacts.py` → `process_all_pdb_contacts.py` →
`get_pdb_parity.py`) was done to check this directly, run standalone
(outside Nextflow, no conda) via a dedicated Python 3.9 venv
(`pip install openbabel-wheel pdbe-arpeggio` — no conda environment
needed, despite `nextflow/envs/arpeggio-env.yaml` being conda-based).

**Structures** (chosen from the plan's own Verification section, EC
confirmed against real SIFTS `pdb_chain_enzyme.tsv.gz` data, not assumed):

| PDB | Protein | EC | Cofactor of interest |
|---|---|---|---|
| 2CPP | Cytochrome P450cam | 1.14.15.1 | heme b |
| 1MTY | Methane monooxygenase hydroxylase | 1.14.13.25 | non-heme diiron center |
| 1JB0 | Photosystem I | 1.97.1.12 | chlorophyll + `[4Fe-4S]` |

**Results, 2 of 3 structures** (1JB0 did not complete — see below, a
confirmed scale limitation, not a bug):

- **2CPP: confirmed positive.** The real bound `HEM` ligand scores
  **0.717 PARITY similarity** (threshold 0.40) against cognate ligand
  `HEA` (heme A) — a row that exists in `cognate_ligands_df.pkl` *only*
  because of the v2 cofactor path (`ligand_source: cofactor`,
  `CofactorDB:22`). Before this work, heme was never a stoichiometric
  reactant/product of the P450cam reaction, so it could not have scored
  against anything. Camphor (the real substrate) separately scores 1.0
  via the pre-existing reaction path, confirming that path is unaffected.
- **1MTY: confirmed limitation, not a bug.** The real bound `FE (III)
  ION` scores 0.0 against every cognate ligand. MMO's actual cofactor is
  a carboxylate-bridged **diiron** cluster; this structure deposits it as
  a bare mononuclear ion, which neither combined source (CoFactor DB's 27
  classes, UniProt's structured entries for this protein) resolves to a
  matching cognate ligand as represented. A genuine, specific gap - not
  general evidence against the mechanism (2CPP's heme match is real
  evidence for it).
- **1JB0: known scale limitation, not a bug — accepted, not pursued
  further.** Re-run to completion and genuinely timed out: PDBe-Arpeggio
  ran for the full 6-hour wall-clock limit on 1JB0's biological assembly
  without finishing, because that assembly is a crystallographic trimer
  (3 symmetric copies of every chain - the selection list contains
  `/A/`, `/A_2/`, `/A_3/`, etc. for each chain), not the monomeric unit -
  a far larger atom count than 2CPP or 1MTY. This is a real limitation of
  running PDBe-Arpeggio standalone on a single machine against a large
  biological assembly, not a defect in the cofactor-coverage work itself.
  Chlorophyll was still expected to fail to match (known gap, see above)
  and the `[4Fe-4S]` cluster was still expected to succeed (confirmed
  present in `cognate_ligands_df.pkl` for this EC via the UniProt path),
  but neither was confirmed against real structure data - this specific
  case remains unconfirmed. Not pursued further this session; a smaller/
  monomeric Photosystem I structure, or a proper Nextflow/cluster
  environment with more time and resources, would be the way to actually
  finish this case if wanted later.

**Two more real bugs found and fixed running this** (same
strict-dtype-coercion/version-drift class as the ones found building
`cognate_ligands_df.pkl` - see the commit history):
`process_pdb_contacts.py` assigning a list into a strictly-string-typed
`xref_db` column, and a `np.where(..., "minor", np.nan)` call that numpy
2.x's stricter dtype promotion now rejects when mixing a string result
with `np.nan` (fixed by using `None` instead).

**This is evidence, not a benchmark.** Two hand-picked structures is
enough to confirm the mechanism works on a real, known case and to
surface one real limitation — it is not a statistically meaningful
measure of match rate across the ~4,843 newly-covered ECs. See the
larger, stratified sample below.

## Larger, stratified PARITY benchmark (2026-10-10)

A proper (if still modest) follow-up to the hand-picked check above:
rather than choosing structures to make a point, this samples the actual
`ligand_source == "cofactor"` population (2,077 rows, 1,833 distinct
ECs — compounds that exist in `cognate_ligands_df.pkl` *only* because of
this work) stratified by cofactor identity, so every distinct cofactor
class (not just the two most common: `[4Fe-4S]` cluster and heme A) gets
at least one real test case. Method:

1. For each of the 1,833 cofactor-only ECs, found its distinct cofactor
   identities (`compound_name`) and, for each (EC, identity) pair, the
   smallest real PDB structure annotated with that EC (SIFTS chain-count
   as a cheap size proxy, to avoid another 1JB0-style stall) via SIFTS
   `pdb_chain_enzyme.tsv.gz` - not every target has a real small
   structure available; 946 of the 1,833 ECs do.
2. Took one EC per distinct cofactor-identity category (random seed 7),
   giving **24 structures across 24 ECs**, covering heme/heme-d variants,
   `[4Fe-4S]`/`[Ni-4Fe-4S]`/`[8Fe-7S]`-type clusters, B12 (as
   methylcobalamin), MGD, SAM, FAD, DPM, PQQ, TPQ/TRQ/LTQ, several
   lanthanide/transition-metal ions, and more.
3. Ran the same standalone pipeline chain as the small check, but with a
   short (15 min, not 6h) per-structure Arpeggio timeout, in a single
   batch rather than per-structure invocations. All 23/24 structures
   cleared Arpeggio (none timed out at this size; one structure, 2nyf,
   had no bound entities in its CIF at all - a data finding, not a
   pipeline failure). One structure (1dd9) produced no domain contacts
   above the cutoff for its ligand and was excluded. `get_pdb_parity.py`
   scored the remaining structures' bound ligands against the real
   `cognate_ligands_df.pkl`.

**Results**: of the 24 targeted ECs, 20 produced usable real structural
data. Among those, 16 (EC, real bound ligand) pairs scored a PARITY
match ≥ 0.40, and **5 of those exist only because of this session's
cofactor-coverage work** (`ligand_source == "cofactor"`, not
`"reaction"` or a mix):

| EC | Real PDB ligand | Score |
|---|---|---|
| 1.1.2.10 | Lanthanum ion | 1.000 |
| 2.5.1.61 | DPM (a dipyrromethane-type cofactor) | 1.000 |
| 1.11.1.10 | Heme | 0.956 |
| 7.1.1.7 | A1JN4 (a heme *d* derivative) | 0.780 |
| 1.2.7.4 | `[Ni-4Fe-4S]` cluster (deposited as `XCC`) | 0.462 |

This is real, diverse evidence beyond the single P450cam/heme case above
- a metal ion, an open-chain porphyrin precursor, two distinct heme
variants, and a nickel-iron-sulfur cluster, each a different EC and a
different real PDB structure, none hand-picked to make this point.

**Honest misses, not glossed over**: `SF4` (a plain iron-sulfur cluster)
at EC 1.7.99.1 scored only 0.333 (below threshold, despite `SF4`-type
clusters being cofactor-covered elsewhere in the dataset - this specific
structure's deposited variant didn't score highly enough). Cobalamin
(B12, as `COB`) at EC 2.8.4.1 scored 0.316 (close, still below). SAM at
the same EC scored only 0.171 despite being a dual-sourced
(`cofactor|reaction`) cognate ligand - the specific deposited structural
variant in that one structure didn't closely match. These are genuine,
structure-specific limitations of PARITY-level matching, not failures of
the EC-coverage mechanism itself (the candidate cognate ligand for each
of these *is* present in `cognate_ligands_df.pkl` - it just didn't score
above threshold against this particular real structure's deposited
form).

**One more real, non-blocking bug found running this**:
`process_all_pdb_contacts.py` crashes later in an unrelated CATH-domain-
interaction-output step when run standalone without `--cache`
(`AssertionError` on a `cath_contacts_mmcif` merge) - but
`bound_entities_to_score.pkl`, the only file `get_pdb_parity.py` actually
needs, is written *before* that step, so it didn't block this result.
Not yet root-caused or fixed - flagged for whoever next touches
`process_all_pdb_contacts.py`'s standalone (non-Nextflow) invocation
path.

**Still not the full picture**: 24 of 1,833 cofactor-only ECs (and only
one structure per EC) is still a small fraction, and the sampling
favoured the smallest available structure per EC rather than a random
one. A genuinely complete benchmark - all cofactor-only ECs, multiple
structures per EC where available - remains open; see the standing
memory to-do.

## Future work

1. **An even larger/fuller PARITY benchmark.** The 24-structure
   stratified sample above meaningfully advances this item (up from 2-3
   hand-picked structures to a real, diverse, non-cherry-picked sample)
   but doesn't close it - it covers 24 of 1,833 cofactor-only ECs, one
   structure each. A fuller version would sample more ECs and more than
   one structure per EC. **Tracked as a standing to-do in memory** (see
   `todo_v2_parity_benchmark` in this project's Claude memory) rather
   than only here, so it survives across sessions.
2. **LLM-assisted extraction for the 742 free-text-only UniProt
   `COFACTOR` rows** — the mechanical fix for chlorophyll specifically
   (see "Known gaps" above).
3. **Root-cause the 9 fully-lost ECs and the ~104 macrocycle-structure
   lost pairs**, if full parity with the pre-v2 dataset is ever required.
4. The 1JB0 (chlorophyll/`[4Fe-4S]`) case specifically — genuinely timed
   out standalone (6h Arpeggio wall-clock limit on its trimeric
   biological assembly), accepted as a known scale limitation and not
   pursued further this session. Revisit with a smaller/monomeric
   Photosystem I structure, or a proper Nextflow/cluster environment,
   if this specific case is wanted later.

See [docs/cofactor_coverage_plan.md](cofactor_coverage_plan.md#remaining-next-steps-not-yet-done)
for the same list with more detail.

## Reproducing this benchmark

```bash
python3 nextflow/bin/benchmark_cognate_ligands.py \
    --old /path/to/old_cognate_ligands_df.pkl \
    --new /path/to/new_cognate_ligands_df.pkl \
    --report_out benchmark_report.txt
```

Re-run this after any future change to the cognate-ligand generation
pipeline (`preprocess_rhea.py`, `preprocess_cofactors.py`,
`get_ec_information.py`) to get a real before/after comparison rather
than assuming the effect of a change.

## See also

- [docs/cofactor_coverage_plan.md](cofactor_coverage_plan.md) — the
  original design/planning document (source analysis, source selection
  rationale, and the session-by-session implementation record this
  summary is drawn from).
- `nextflow/bin/preprocess_cofactors.py`,
  `nextflow/bin/benchmark_cognate_ligands.py` — implementation and
  benchmark script.
