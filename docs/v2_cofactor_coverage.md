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
  — chlorophyll/heme-biosynthesis-pathway intermediates specifically,
  plausibly an RDKit sanitization edge case on unusual valence states
  (characterised, not yet root-caused).
- 270 are ordinary, otherwise-well-resolved compounds missing only for
  specific ECs — consistent with normal upstream Rhea/KEGG reaction-
  equation revision over the ~2-year gap (confirmed example: NADH is
  missing specifically for EC 1.1.1.96, but still resolves correctly for
  630 other ECs — the reaction equation for that one EC now cites NAD
  where it previously cited NADH).
- 9 ECs lost entirely (`1.14.14.140, 2.1.1.86, 2.4.1.129, 2.7.11.27,
  3.1.8.2, 4.2.1.78, 4.3.3.2, 4.3.3.3, 4.3.3.4`) — not yet individually
  root-caused; small enough (9 of 6,214) not to block on, but worth
  revisiting if full v1 parity is ever required.

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
