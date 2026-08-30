# Plan: Closing the Cofactor Coverage Gap

## The problem

Cognate ligands in ProCogGraph are currently derived exclusively from
**reaction equations**: `enzyme.dat` EC entry → KEGG enzyme record → KEGG
reaction → reactant/product compound codes → ChEBI/PubChem SMILES
(`nextflow/bin/get_ec_information.py:415-589`, cross-checked against Rhea
in `nextflow/bin/preprocess_rhea.py`). A compound only becomes a candidate
cognate ligand if it is explicitly written as a reactant or product in the
balanced equation for that EC number.

This is fine for substrates and products, and for cofactors that *are*
stoichiometric participants (e.g. `NAD+ + substrate ⇌ NADH + product` —
NAD is written into the equation, so it's captured). It structurally fails
for cofactors that participate **catalytically rather than
stoichiometrically** — prosthetic groups, electron/light carriers,
structural metal ions — because these are never written as a reactant or
product of the EC's net reaction at all.

Your thesis identifies exactly this (p.172): Chlorophyll A is the single
most frequently unmatched ligand in ProCogGraph, "due to its lack of
annotation as a cofactor in reaction schemes" — chlorophyll is never
consumed or produced by a photosystem's EC reaction, so no amount of
reaction-database improvement will surface it via the current pipeline
path.

`enzyme.dat` itself doesn't help either: ExPASy dropped the structured
`CF` cofactor field years ago; cofactor mentions today only exist as
unstructured free text inside `CC` comment lines (`utils.process_ec_records`,
`nextflow/bin/utils.py:17-37`, currently only extracts `ID`/`DE`/derived
`TRANSFER`) — not reliably machine-parseable.

### An existing, separate mechanism that this plan is not about

The pipeline already has a cofactor-related step, at
`get_ec_information.py:699-719`: after `cognate_ligands_df` is built from
the reaction-equation path, it joins each ligand's ChEBI ID against
`chebi_relations.tsv`'s `has_role` relations for four ChEBI role classes
(`CHEBI:23357` cofactor, `23354` coenzyme, `26348` prosthetic group,
`26672` siderophore) and stamps an `isCofactor` column. This **labels**
ligands that are already in `cognate_ligands_df` — it does nothing to add
ligands that never entered the table because they're not reaction
participants. It's a labelling mechanism for the existing set, not a
coverage extension. This plan is about the latter; the two are
independent and should stay that way, but the naming overlap (`isCofactor`
vs. the new `ligand_source` value below) is worth keeping straight when
implementing.

(Also checked: a `parse_brenda_json_generic_reaction` function exists at
`get_ec_information.py:323` but is dead code — never called anywhere.
BRENDA is not currently wired into this pipeline despite that function's
presence.)

## Sources considered

### Ruled out: PDBe RelLig bulk TSVs (original plan's primary source)

The original version of this plan proposed PDBe's **RelLig**
(`https://github.com/PDBeurope/rellig`) bulk PDBeChem v2 output
(`interacting_chains_with_ligand_functions.tsv`,
`pdb_bound_molecules.tsv`) as the primary source, on the assumption these
files carried an EC number per row. **Verified against the real file and
this is wrong** — pulled the actual header + rows:

```
PDBID  Chain_Symmetry  BestUnpAccession  LigandID  bmID  LigandType  inchikey
ProteinNameUniprot  ProteinNamePdb  OrganismScientificNameUniprot  ...  annotation
101m   A  P02185  HEM  bm2  CCD  KABFMIBPWCXCRK...  Myoglobin  Myoglobin  ...  reactant-like
```

No EC column. RelLig's actual EC-bearing output is a **separate per-ligand
JSON format** (`<ligand_id>_cofactor_annotation.json`), produced by
*running* `pdberellig cofactors --cif ...` per CCD component — not a
static bulk download. Chasing this further would mean either running the
pipeline ourselves against the full CCD, or finding pre-generated JSON
output that may not exist as a bulk artifact at all. Not pursued further.

### What RelLig's source code did turn up, and is being reused instead

Investigated `pdberellig`'s own source
(`pdberellig/core/cofactors.py`, `pdberellig/data/cofactors/`) directly.
It ships two small, static, redistributable (Apache-2.0) data files that
are exactly the EC↔cofactor mapping this plan actually needs, with none of
the bulk-TSV baggage:

- **`cofactor_ec.csv`** — 3,915 rows, `EC_NO → COFACTOR_ID`, sourced from
  two curated inputs: `cofactor_db_2010` (the **CoFactor database**,
  Fischer/Holliday/Thornton, *Bioinformatics* 2010 — a manually curated
  catalogue of organic cofactors and the EC numbers known to use them) and
  `brenda` (supplementary associations from BRENDA — note this is RelLig's
  own upstream use of BRENDA, still not something this codebase depends on
  directly).
- **`cofactors_details.json`** — 27 distinct cofactor classes (IDs 1–28,
  one retired), each with a representative PDB CCD code (id 4 → `NAD`, id
  22 → `HEA`, id 7 → `PLP`, id 11 → `B12`, etc. — full list: TPP, FAD,
  FMN, NAD, pantetheine, CoA, PLP, glutathione, biotin, folate, B12,
  ascorbate, menaquinone, ubiquinone, molybdopterin, tetrahydrobiopterin,
  a mycofactocin-type cofactor, SAM, coenzyme F430, coenzyme M, heme A,
  deazaflavin, PQQ, TPQ, TRQ, lipoic acid).

Since `ccd_cif` is already a pipeline input, resolving a representative
CCD code to a SMILES needs no new source at all.

**Known gap in this source, confirmed and not fixable by construction**:
chlorophyll is not among the 27 classes. CoFactor DB (2010) is a
relatively small, central-metabolism-focused curated set — it doesn't
cover chlorophyll, and doesn't cover bare structural metal ions beyond
what's implicit in a couple of its classes. This source alone does **not**
solve the thesis's headline motivating example.

### Added: UniProt `COFACTOR` annotation, structured entries

Originally scoped in this plan as a secondary/fallback source (per-accession,
chunked REST calls). Re-investigated as a first-class combined source
instead, and it's better suited than that framing suggested:

- **Bulk-downloadable**, not per-accession: confirmed the `/uniprotkb/search`
  (or `/uniprotkb/stream` for larger pulls) endpoint supports
  `fields=accession,ec,cc_cofactor&format=tsv` directly —
  `reviewed:true AND cc_cofactor:* AND ec:*` returns 111,819 rows in one
  bulk pull, no per-accession chunking needed.
- **99.3% of those rows (111,077) already have a machine-parseable
  `Name=...; Xref=ChEBI:CHEBI:NNNN;` structure** — confirmed by regex
  extraction against the real bulk TSV. Only 742 rows (0.7%) are
  free-text-only (a bare `Note=` with no `Xref`).
- **Chlorophyll's specific entries fall into that 0.7% free-text bucket** —
  confirmed directly: `cc_cofactor:"ChEBI:CHEBI:18230"` (chlorophyll a's
  real ChEBI ID) returns **zero** reviewed UniProt entries. The actual
  annotation for e.g. Photosystem I (`P56766`, EC 1.97.1.12) is:
  ```
  COFACTOR: Note=P700 is a chlorophyll a/chlorophyll a' dimer, A0 is one
  or more chlorophyll a, A1 is one or both phylloquinones and FX is a
  shared 4Fe-4S iron-sulfur center.
  ```
  No `Name=`/`Xref=` — not mechanically resolvable to a SMILES. So even
  the broader UniProt source doesn't mechanically solve chlorophyll; it's
  annotated qualitatively (variable stoichiometry, mixed a/a' forms), not
  as a clean 1:1 cofactor identity, in both curated sources checked.

## Why combine both rather than pick one

Measured real overlap between the two sources (terminal-EC level, exact
string match, before any broadcast):

| Source | Distinct ECs covered |
|---|---|
| CoFactor DB 2010 + BRENDA (`cofactor_ec.csv`) | 2,760 |
| UniProt `COFACTOR`, structured only | 3,284 |
| **Overlap** | 1,201 |
| **CoFactor DB-only** (UniProt has no cofactor annotation at all here) | 1,559 |
| **UniProt-only** (outside CoFactor DB's 27-class scope entirely) | 2,083 |
| **Combined unique** | **4,843** |

Confirmed this is real complementary coverage, not redundancy, with a
concrete example: **EC 1.1.1.10** (D-xylulose reductase, a textbook
NADP-dependent enzyme) is in CoFactor DB's list, but every reviewed
UniProt entry for it has a **completely empty** `Cofactor` field (checked
live: `Q7Z4W1`, `Q91X52`, `Q21929` all blank). CoFactor DB/BRENDA encode
EC-class-level literature knowledge that individual UniProt curators never
entered per-accession; UniProt in turn catches specific, curator-verified
cofactor identities (particularly metal centers — 109 distinct ChEBI
identities vs. CoFactor DB's 27 classes) that fall outside CoFactor DB's
fixed, 2010-dated scope. Use both, unioned.

## The broadcast problem, and its fix

Not all of UniProt's EC values are fully resolved to 4 digits — 196 of the
3,284 (6.0%) are partial (`N.N.N.-`, `N.N.-.-`, or `N.-.-.-`). A naive
exact-string join (the same pattern `get_pdb_parity.py` already uses for
reaction-derived ligands) silently drops all 196, undercounting coverage:
exact-match-only gives just **62.0%** of the pipeline's 6,753 real terminal
ECs (4,184/6,753).

**First attempt at "broadcast the wildcard down to matching terminal ECs"
produced a bogus 100% coverage number** — traced this to the 7 class-level
wildcards (`1.-.-.-` through `7.-.-.-`) each matching literally every
terminal EC in that class (e.g. `1.-.-.-` → all 6,753... no, all of class
1). Broadcasting at that granularity is chemically wrong — it would claim
catalase (`1.11.1.6`, heme-dependent) shares a cofactor with an unrelated
NAD-dependent dehydrogenase just because both are oxidoreductases.
Breaking down by granularity confirmed the danger is real:

| Partial-EC granularity | Count | Terminal ECs it would broadcast to |
|---|---|---|
| `N.-.-.-` (class-level) | 7 | 6,753 — literally everything |
| `N.N.-.-` (subclass-level) | 33 | 5,578 |
| `N.N.N.-` (subsubclass-level) | 156 | 6,193 |

**Rule adopted**: only broadcast **subsubclass-level (`N.N.N.-`)** partial
ECs — enzymes sharing all three EC digits genuinely tend to share cofactor
chemistry (e.g. `1.1.1.-` is the NAD(P)-dependent CH-OH oxidoreductase
subsubclass). **Class- and subclass-level partial ECs (40 of the 196) are
dropped from the cofactor table entirely** — there's no chemically
defensible way to narrow "some oxidoreductase, unknown subclass" down to
specific terminal children.

With that rule, combined coverage of the pipeline's real terminal EC list:

- Exact match only (CoFactor DB + UniProt 4-level): **62.0%** (4,184/6,753)
- \+ safe subsubclass-level broadcast: **96.2%** (6,496/6,753)

**Where the broadcast happens**: at `cognate_ligands_df` *build* time, not
at match time. For each UniProt row that only resolved to `N.N.N.-`, look
up the pipeline's own terminal EC list (`ec_records_df.TRANSFER.unique()`,
already computed in `get_ec_information.py`'s `main()`) and expand that
one row into one row per real terminal EC sharing that subsubclass prefix,
each carrying the same cofactor ChEBI ID/SMILES — identical shape to every
existing reaction-derived row. This means `cognate_ligands_df.entry` never
contains a wildcard, so `get_pdb_parity.py:121`'s existing
`cognate_ligands_df.entry.isin(ec)` join needs **zero modification** — it
simply sees more rows for more ECs. No new match-time logic anywhere
downstream. This was the whole point of preferring this approach over the
RelLig structure-instance path: cofactor rows are ordinary
`cognate_ligands_df` rows, keyed by real EC, from the start.

## Known remaining gap: chlorophyll and other free-text-only cases

742 UniProt rows (0.7% of the reviewed cofactor+EC set) have only a
free-text `Note=` with no structured `Xref=ChEBI:...` — chlorophyll's
photosystem entries are in this bucket. Not mechanically resolvable to a
SMILES from either source as currently curated.

**Deferred, explicitly agreed**: revisit this bucket with LLM-assisted
extraction (read the free-text note, propose a ChEBI ID / SMILES, flag for
review) as a **separate follow-on pass** once the main combined table
(CoFactor DB + UniProt structured, with subsubclass broadcast) is built
and integrated. Don't block the main implementation on this — it's a
small, well-bounded tail cleanup (742 rows), not core to getting the
combined source working.

## Implementation steps

1. **Vendor `cofactor_ec.csv` and `cofactors_details.json`** from
   `pdberellig` (Apache-2.0, redistribution permitted) into the reference
   data set, alongside a script step resolving each `COFACTOR_ID`'s
   representative CCD code to canonical SMILES via the existing
   `ccd_cif`-parsing path.
2. **Bulk-pull UniProt structured cofactor data**: add
   `https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+cc_cofactor:*+AND+ec:*&fields=accession,ec,cc_cofactor&format=tsv`
   (or the equivalent `search` pagination) to the reference data manifest
   pattern already established in `docs/reference_data_download_plan.md`.
   Parse `Name=...; Xref=ChEBI:CHEBI:(\d+)` per row; keep only rows with a
   match. Resolve ChEBI ID → SMILES via the same ChEBI resolution path
   `get_ec_information.py` already uses for reaction-derived ligands.
3. **Split UniProt rows by EC completeness**: exact 4-digit ECs pass
   through directly; `N.N.N.-` rows get expanded against the pipeline's
   own terminal EC list (`ec_records_df.TRANSFER.unique()`) at build time;
   `N.-.-.-` and `N.N.-.-` rows are dropped (logged, not silently
   discarded, so the drop is auditable).
4. **Union CoFactor DB rows + expanded UniProt rows**, tag both with a new
   `ligand_source = "cofactor"` value (existing reaction-derived rows get
   `ligand_source = "reaction"` for consistency — this column doesn't
   exist yet and needs adding to the reaction-derived rows too, not just
   the new ones). Concat into `cognate_ligands_df` alongside the existing
   Rhea/KEGG/ChEBI/PubChem/GlyTouCan frames at
   `get_ec_information.py:675`.
5. **Do not modify `get_pdb_parity.py` or `produce_neo4j_files.py`** — the
   existing EC-join and PARITY-scoring logic should work unchanged against
   the enlarged `cognate_ligands_df`, since cofactor rows are
   indistinguishable in shape from reaction rows except for the new
   `ligand_source` tag.
6. **(Follow-on, separate piece of work)**: LLM-assisted extraction for
   the 742 free-text-only UniProt rows, chlorophyll included.

## Verification

- After implementation, re-check the specific case cited in your thesis
  — note per the gap analysis above, chlorophyll is **expected to still be
  missing** after step 1–5 alone; don't treat its absence as a bug until
  the follow-on LLM-extraction pass (step 6) is done.
- Confirm `cognate_ligands_df.entry` contains zero wildcard/partial EC
  values after the broadcast step — this is a hard invariant the
  downstream `isin()` join depends on.
- Spot-check a handful of well-characterised cofactor-dependent enzymes
  outside chlorophyll (e.g. cytochromes with heme, the *Methylococcus
  capsulatus* MMO example already discussed in your thesis's future-work
  section, p.175) to confirm the new cofactor-origin ligands match known
  biology.
- Confirm reaction-origin cognate ligands are completely unaffected (this
  is a strict addition) — existing PARITY scores and cognate ligand counts
  for already-mapped structures should be unchanged, and the new
  `ligand_source` column should backfill to `"reaction"` for every
  pre-existing row.
- Compare before/after counts of "most frequently unmatched ligands" (the
  analysis behind Table 3.3 / the Discussion's unmatched-ligand
  discussion) for the ~4,843 newly-EC-covered enzymes.

## Benchmark results

Real end-to-end build (2026-08-30), benchmarked against the last
pre-cofactor-coverage `cognate_ligands_df.pkl` (2024-07): **+87%
EC-ligand pairs (41,785 → 78,163), +1,102 net ECs covered (6,214 →
7,316)**, with the "strict addition" design intent confirmed (40,282
pairs unchanged) and a small, fully-categorised set of differences from
the old build (3.6%, mostly RDKit macrocycle-sanitization edge cases and
normal upstream Rhea/KEGG data drift, not a regression in this plan's own
logic).

Full methodology, the complete results table, the lost-pair
categorisation, and reproduction instructions now live in
**[docs/v2_cofactor_coverage.md](v2_cofactor_coverage.md)** — written as
the first-class, citable record of this work (for the v2 docs/paper),
rather than duplicated here where it would drift out of sync.

## Status

**Implemented and benchmarked** — see
[docs/v2_cofactor_coverage.md](v2_cofactor_coverage.md) for the full
writeup. What was originally scoped as "suggested order of work" is done:
vendoring/parsing the CoFactor DB tables, the UniProt bulk pull and
exact/broadcast/drop EC-completeness split, the `ligand_source` column
and union into `cognate_ligands_df`, and a real before/after comparison
against the last pre-cofactor-coverage build.

## Remaining next steps (not yet done)

1. **Real, larger-scale PARITY-based validation against PDB structures.**
   A *small* version of this is now done — see
   [docs/v2_cofactor_coverage.md](v2_cofactor_coverage.md#small-scale-parity-validation-2026-08-3031):
   `get_pdb_parity.py` run standalone against 2 of 3 hand-picked
   structures (P450cam/heme, MMO/diiron; Photosystem I incomplete),
   confirming a real cofactor-origin PARITY match on P450cam's heme
   (0.717, clears threshold) and a real, specific gap on MMO's bare Fe³⁺
   ion (0.0, not matched). That's evidence the mechanism works, not a
   statistically meaningful measure across the ~4,843 newly-covered ECs —
   a proper sampled benchmark is still open. **Tracked as a standing
   to-do in memory** (`todo_v2_parity_benchmark`), so it isn't lost
   between sessions.
2. **LLM-assisted extraction for the 742 free-text-only UniProt
   `COFACTOR` rows** (chlorophyll's actual fix - see "Known remaining
   gap" above). Explicitly deferred, not started.
3. **Investigate the 9 ECs and ~104 macrocycle-structure pairs lost
   relative to the pre-cofactor-coverage baseline**, if full parity with
   the old dataset is ever required (see
   [docs/v2_cofactor_coverage.md](v2_cofactor_coverage.md)'s "Known
   gaps" for the specific EC list and categorisation). Not blocking -
   small (9 of 6,214 ECs) - but not yet root-caused either.
4. Finish the stopped 1JB0 (Photosystem I) run — the chlorophyll/`[4Fe-4S]`
   case — likely as part of whichever of the above happens first.
