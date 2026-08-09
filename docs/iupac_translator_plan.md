# Plan: IUPAC-Extended → IUPAC-Condensed Translator (for offline glycan→SMILES)

## Context (read this first if you're picking this up cold)

ProCogGraph's cognate-ligand generation pipeline resolves glycan
structures to SMILES via a chain of live external APIs (GlyTouCan,
GlycoSmos, CSDB). Investigation earlier in this project (see
`docs/cognate_ligands_generation_plan.md`, Part 1 "Question 5" and
"Question 6") found:

- One of those live dependencies (GlyTouCan's `gtcid2seqs` API) is
  **currently broken in production** — it stopped returning the field the
  pipeline code depends on (`glycoct`), confirmed via a 100-sample test:
  0% success live, vs. ~54% success using an offline library (`GlyLES`)
  instead.
- `GlyLES` (`pip install glyles`, MIT license) converts **IUPAC-condensed**
  glycan notation (e.g. `GlcNAc(b1-4)GlcNAc`) to SMILES, entirely offline,
  no API calls.
- `GlyLES` cannot take **WURCS** as input via its public API — but WURCS
  is what's actually available in two places that matter: (a) GlyTouCan's
  API still reliably returns `wurcs` even though `glycoct` is broken, and
  (b) PDB structures themselves only ever carry WURCS/GLYCAM/LINUCS for
  their bound carbohydrate entities
  (`_pdbx_entity_branch_descriptor.type`), never IUPAC-condensed —
  confirmed against the mmCIF dictionary.
  - **Checked and ruled out in this session**: GlyLES actually ships its
    own WURCS grammar/lexer/parser internally (`glyles/wurcs/`,
    `WURCSGlycan` class, exported at package top-level) — this looked like
    it might be an undocumented shortcut that bypasses the whole
    translator problem. It doesn't work: the lexer chokes on standard
    WURCS syntax (`extraneous input 'a'...`), and even if parsing
    succeeded, the actual tree-walk/merge-to-SMILES step is **commented
    out in the source** (`wurcs_glycan.py`, the `self.parse_tree = ...` /
    `Merger(...).merge(...)` lines are dead code behind a comment). It's
    an unfinished upstream feature, not usable today. This confirms the
    glypy→translator→GlyLES route below is still the right path, not a
    detour.
- `glypy` (`pip install glypy`, Apache-2.0, actively maintained,
  `mobiusklein/glypy` on GitHub) **can** parse WURCS offline and export
  various formats including IUPAC — verified its GlycoCT export is
  byte-identical to the live GlycoSmos API's output for a real test case.
  But glypy's IUPAC *text* output is a **different dialect**
  ("IUPAC-extended", e.g. `b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc`) than what
  GlyLES needs ("IUPAC-condensed", e.g. `GlcNAc(b1-4)GlcNAc`). Feeding
  glypy's text output directly into GlyLES fails to parse.
  - **Superseded finding, kept for history**: the first prototype (v1,
    see "Superseded: v1" section below) worked by regexing glypy's *text*
    IUPAC-extended output. It proved the algorithm but was explicitly
    flagged as not the final approach. **This session replaced it with
    v2**, which walks glypy's *structured object model* instead of its
    text output — see "What's now built: v2" below. glypy's text IUPAC
    export is no longer used by the translator at all.
- The other realistic offline alternative, `GlycanFormatConverter` (which
  does WURCS→IUPAC-condensed directly), is a **dead end** — its build
  depends on an artifact (`org.glycomedb:residuetranslator:0.1-b12`) that
  has been lost from its only host server (confirmed 404, not just an
  access/cert problem; checked Maven Central, GitHub, and the Wayback
  Machine — no copy exists anywhere).

**This translator is the missing piece.** If it works reliably, the full
chain `WURCS → glypy (offline) → [this translator] → GlyLES (offline) →
SMILES` becomes possible with **zero external API calls**, fixing the
currently-broken cognate-ligand glycan path and giving a fully offline
alternative for PDB bound-entity glycan resolution too.

## What's now built: v2 (structured object-model translator, this session)

**Status: working prototype, validated against 4 real cases including all
previously-proven ones plus new residue types v1 never covered. Not yet
integrated into the pipeline, not yet scale-tested.**

### Why v2 replaces v1's approach entirely

v1 (see "Superseded: v1" below) worked by regexing glypy's **text**
IUPAC-extended output and had a 5-entry hand-built `RESIDUE_MAP`. The
plan's own "what's not done" list flagged two things to check before
expanding that table by hand: whether GlyLES exposes its own vocabulary
programmatically, and whether glypy exposes a structured object model
instead of just text. Both turned out to be true, and **using them
directly is a much smaller and more robust implementation than a flat
lookup table**, for a specific reason found while investigating GlyLES's
own grammar (`glyles/iupac/IUPAC.g4`):

**GlyLES's IUPAC-condensed input is not a flat list of literal residue
names.** Its grammar decomposes each residue as `modi* saci+ modi*` —
a root sugar token (`saci`, one of 34 literal names: `Glc`, `Man`, `Gal`,
`Neu`, `Kdo`, `Fuc`, ... — matches glypy's own base-sugar identities
almost 1:1) plus zero or more `modi` modifier tokens, each of which is an
optional position number + an optional bridging-atom letter (`C`/`N`/`O`/
`P` — the atom the modification attaches through) + a functional-group
code (`Ac`, `Gc`, `S`, `Me`, `P`, ... — from a fixed, enumerable list in
the grammar). So `GlcNAc` is not an atomic string GlyLES recognizes by
name — it's literally `Glc` + modifier(`bridge=N`, `fg=Ac`, implied
position 2), and GlyLES accepts the fully-explicit form `Glc2NAc` as
identical (verified: byte-identical SMILES for both). Likewise `Neu5Ac`
decomposes structurally into the base sugar `Kdn` + modifier(`5`, `N`,
`Ac`) — verified: `Kdn5NAc`, `Neu5Ac`, and `NeuAc` all produce
**byte-identical SMILES** through real GlyLES. This exactly mirrors how
glypy itself represents these residues internally (see below) — the two
libraries' internal models line up almost exactly once you go past their
respective text formats.

This means the translation problem is not "map every named glycan residue
(dozens, combinatorially exploding with modifications) to its GlyLES
string" — it's "map glypy's ~15 distinct base-sugar structures and ~35
distinct substituent types (both fixed, enumerable, already defined in
each library's own source) to GlyLES's grammar tokens." Much smaller
surface area, and correct by construction for combinations neither
library's authors thought to name explicitly (e.g. a hypothetical
disulfated, diacetylated hexose — v1's flat-table approach could never
have covered that; v2 does, automatically, because it composes rather
than looks up).

### How v2 works

1. **Base sugar identity**: for each glypy `Monosaccharide` node, clone it,
   strip its substituents (`node.substituents()` → drop each), and call
   `glypy.io.nomenclature.identity.identify()` on the stripped clone. This
   matches the residue's ring/stereochemistry/backbone against glypy's own
   ~100-entry named-residue registry (`glypy.monosaccharides`) and returns
   a name — generically, for *any* residue in that registry, not just ones
   someone explicitly added to a lookup table. Verified this correctly
   reduces `GlcNAc` → `Glc`, `NeuAc`/`NeuGc` → `Kdn` (their true structural
   base, per glypy's own model), and directly returns the right literal
   shortcut names for deoxyhexoses and acids that GlyLES also recognizes
   as single tokens (`Fuc`, `Rha`, `Qui`, `GlcA`, `IdoA`, `Kdo` all come
   back correctly with zero special-casing needed — they're already
   distinct entries in glypy's registry with the right stereochemistry
   baked in).
   - **One real quirk found and handled**: `identify()` is
     anomer-sensitive (it matches the query's anomer against each
     registry candidate's own fixed anomer), and glypy's registry has
     exactly one entry, `aMan`, whose *name* bakes in a specific anomer
     rather than being anomer-generic — so a real alpha-configured
     mannose residue matches `aMan`, not `Man`. Since we already track
     anomer separately (from `node.anomer`, used at the linkage-emission
     step), we strip this back off via a 1-entry alias table
     (`ANOMER_PREFIXED_ALIASES`). Confirmed by inspection that `aMan` is
     the only such entry in the whole ~100-name registry — not expected
     to need extending, but documented in the code in case a future glypy
     version adds more.
2. **Substituents**: `node.substituents()` returns structured
   `(position, Substituent)` pairs. Each `Substituent.name` (e.g.
   `n_acetyl`, `sulfate`, `phosphate`, `methyl`) is looked up in a small
   static table (`SUBSTITUENT_MAP`) mapping it to a `(bridge_letter,
   FG_code)` pair, then formatted as `{position}{bridge}{FG}` (e.g.
   position 2, bridge `N`, FG `Ac` → `"2NAc"`). glypy's own substituent
   vocabulary is fixed and small (~35 entries total, enumerated in
   `glypy/structure/substituent.py`'s two `DerivatizePathway` dicts) — so
   this table only needs to grow to that bound, not to cover every
   modified-residue *name* combinatorially.

   **Important correction to how this table was built (found and fixed
   this session, worth flagging so it isn't repeated)**: the first pass
   at this table was built by *reasoning about the grammar* (`bridge* fgi`
   → assume any bridge+FG combination composes correctly) rather than by
   testing every entry against real GlyLES output. That assumption is
   **false**. Inspecting GlyLES's actual reactor code
   (`glyles/glycans/mono/reactor.py`, the `functional_groups` dict) shows
   each FG code maps to a **hand-written SMILES fragment that already
   includes its own attachment atom** — e.g. `"Me": "OC"` (methyl attaches
   via its own oxygen, for an ether), `"Ac": "C(=O)C"` (acetyl attaches
   directly via its own carbonyl carbon, no oxygen). Composing a bridge
   atom in front of these is only chemically correct when the FG
   fragment's own leading atom is compatible — e.g. `N` + `Ac`'s
   `"C(=O)C"` → `NC(=O)C`, a correct amide. But `N` + `Me`'s `"OC"` →
   `NOC`, which is an **N-methoxyamine, not an N-methylamino group** —
   wrong chemistry, found by literally testing `Glc2NMe` and reading the
   output SMILES rather than trusting the grammar shape. **Every entry
   below was verified this way** (built a minimal test glycan bearing the
   substituent, ran it through real GlyLES, and read the resulting SMILES
   to confirm it matches the intended chemistry) — none are guesses.

   - **Mapped and verified correct** (real GlyLES output checked bond-by-bond):
     `n_acetyl` (`N`,`Ac`), `n_glycolyl` (`N`,`Gc`), `n_sulfate` (`N`,`S`),
     `n_amidino` (`N`,`Am`), `sulfate` (`None`,`S`), `phosphate`
     (`None`,`P`), `methyl` (`None`,`Me`), `acetyl` (`None`,`Ac`), `amino`
     (`N`,`None` — a bare bridge atom is itself a valid grammar terminal,
     `fgi: COUNT | bridge | FG`, confirmed via `Glc2N` → correct free
     amine), `glycolyl` (`None`,`Gc` — the non-N/ester form, distinct from
     `n_glycolyl`), `formyl` (`None`,`Fo` — ditto, ester form only),
     `ethyl` (`None`,`Et`), `ethanolamine` (`None`,`Etn`),
     `phospho_ethanolamine` (`P`,`Etn`), `phospho_choline` (`P`,`Cho`),
     `succinate` (`None`,`Suc`), `bromo` (`None`,`Br`), `chloro`
     (`None`,`Cl`), `fluoro` (`None`,`F`), `iodo` (`None`,`I`).
   - **Confirmed unsupported by GlyLES itself, not just unmapped** (tested
     and found to produce wrong chemistry, or ruled out structurally —
     don't waste time re-deriving these): `n_methyl` and `n_formyl` (the
     naive N-bridge composition gives N-methoxyamine / N-formyloxyamine,
     both wrong; `reactor.py` has one hardcoded override,
     `"NFo": "NC=O"`, for a correct N-formyl fragment, but it wasn't
     reached via the straightforward `{pos}N Fo` grammar input tested
     here — worth a deeper look at `reactor.py`'s `set_fg`/name-matching
     logic if N-formyl coverage ever actually matters for real data);
     `n_succinate` (same wrong-composition pattern as n_methyl/n_formyl);
     `thio` (no sulfur bridge token exists in the grammar at all —
     `bridge: CARBON | NITROGEN | OXYGEN | PHOSPHOR`, so this is a real
     grammar limitation, not a documentation gap); `pyrophosphate` and
     `triphosphate` (no corresponding entries exist anywhere in
     `reactor.py`'s `functional_groups` dict); `pyruvate` (GlyLES's `Pyr`
     FG code only encodes a simple single-position ester, verified via
     `Glc6Pyr`; real biological pyruvylation is a cyclic acetal bridging
     *two* ring positions — a two-position `Glc4,6Pyr`-style input was
     tried and GlyLES returned empty SMILES, i.e. it doesn't support the
     real substituent this glypy name refers to); `imino` and
     `dimethylamine`/`n_dimethyl` (no working single-token encoding
     found).
   - **Not attempted, lower priority** (rare in practice, or need
     dedicated non-substituent-table handling rather than a
     `SUBSTITUENT_MAP` entry): `anhydro` (glyles' grammar has a wholly
     separate `modi` alternative for this, `NUM DASH ANHYDRO DASH`, a
     literal keyword-based rule, not a bridge+FG composition — would need
     its own code path, not a table row), `lactone`,
     `diphospho_ethanolamine`, `hydroxymethyl` (would need glyles'
     separate poly-carbon/`carb` branch grammar, meant for long
     branched-chain sugars, not a simple bridge+FG modifier).
3. **Tree structure**: glypy's `Glycan` object exposes real parent/child
   `Link` objects per node (`node.links`, each with `.parent`, `.child`,
   `.parent_position`, `.child_position`) — no text-bracket matching
   needed. Walked recursively from `glycan.root` (the reducing end).
   GlyLES's own bracket convention (confirmed against real GlyLES output
   for the Man3GlcNAc2 case) is: at a branch point, the *last*-positioned
   child continues the main chain unbracketed, and every other child is
   wrapped in `[...]` and emitted before it. v2 replicates this by
   sorting a node's children-links by parent position and treating the
   last one as the main continuation.
   - **Verified beyond the 2-branch case this session**: built a
     synthetic 3-branch node (one parent, three children) and confirmed
     v2's stacked-bracket output (`Man(a1-4)[Fuc(a1-2)][Xyl(b1-3)]Gal`)
     parses through real GlyLES to valid, non-empty SMILES. (An earlier
     attempt at this test produced empty SMILES from GlyLES — turned out
     to be a test-setup mistake, an unset/`?` anomer on one synthetic
     residue, not a translator bug; fixing the test's anomers fixed it.)

### Validated end-to-end (WURCS/structure → glypy → v2 translator → GlyLES → SMILES)

1. **Chitobiose** (same case v1 proved):
   `WURCS=2.0/1,2,1/[a2122h-1b_1-5_2*NCC/3=O]/1-1/a4-b1` →
   `Glc2NAc(b1-4)Glc2NAc` → SMILES **byte-identical** to v1's validated
   output for `GlcNAc(b1-4)GlcNAc`.
2. **Man3GlcNAc2 N-glycan core** (same case v1 proved, the branched one):
   → `Man(a1-6)[Man(a1-3)]Man(a1-4)Glc2NAc(b1-4)Glc2NAc` → SMILES
   **byte-identical** to v1's validated output.
3. **Neu5Ac(a2-3)Gal** (new — v1 never covered sialic acids; v1's table
   only had 5 entries and none were acidic/nonulosonic sugars): built
   directly from glypy's named-residue objects (no real WURCS string
   needed for this spot-check) → `Kdn5NAc(a2-3)Gal` → valid SMILES,
   confirmed structurally equivalent to feeding GlyLES `Neu5Ac(a2-3)Gal`
   directly (byte-identical SMILES both ways).
4. **Synthetic 3-branch case** (Fuc/Xyl/Man off one Gal) — confirmed the
   bracket-stacking logic generalizes past 2 branches, see above.

### Code — now a real module: `nextflow/bin/wurcs_to_iupac.py`

Moved out of the scratchpad prototype into an actual module in the repo
(not yet committed — still an untracked working-tree file, see
`docs/SESSION_HANDOFF.md`'s repo-state section). It exposes a single
public function, `translate(wurcs_string) -> str`, and carries the
benchmark comparison from §2/§2b/§2c in its module docstring so the
numbers travel with the code. **Not yet wired into the pipeline** —
`get_ec_information.py` and `process_all_pdb_contacts.py` don't import it
yet (that's still §3 below), and `glypy`/`glyles` haven't been added to
`nextflow/envs/procoggraph.yaml` yet either.

The module's structure, in brief (full implementation and all the
rationale comments for each non-obvious choice live in the module file
itself, not duplicated here to avoid the two drifting):

- `translate(wurcs_string) -> str` — the one public entry point.
- `_base_sugar_name(node)` — identifies a residue's base sugar generically
  via `glypy.io.nomenclature.identity.identify()` on a substituent-stripped
  clone, with an anomer-override retry loop (handles undefined-anomer
  residues, e.g. free reducing ends) and a furanose-ring `f`-suffix rule.
- `SUBSTITUENT_MAP` — ~19 entries, glypy substituent name → GlyLES
  bridge+FG tokens, each individually verified against real GlyLES output.
- `_emit`/`_emit_edge` — walks glypy's real parent/child `Link` graph to
  build GlyLES's bracket-and-linkage text, handling glypy's
  `UnknownPosition` sentinel (`-1` → GlyLES's `?`) and same-position
  sibling ties that a naive sort would crash on.

A real operational gotcha carried forward from v1's testing still applies:
feeding GlyLES a malformed/half-translated string can hang it for 30+
seconds instead of failing cleanly (suspected ANTLR pathological
backtracking). Any production use of GlyLES needs a subprocess timeout,
regardless of which translator version feeds it. The scale test below
was run with a 15-second `SIGALRM` wrapper around every `glyles.convert`
call for exactly this reason - no hangs were hit in practice on this
sample, but don't skip the safeguard in production.

## What's NOT done — the actual remaining work

### 1. Substituent vocabulary completeness — mostly done this session

`SUBSTITUENT_MAP` now has 19 entries, **all individually verified against
real GlyLES output** (not grammar-shape guesses — see the correction
noted above, where the first-pass guessed table had two wrong entries
caught by actually testing). Of glypy's ~35-entry total substituent
vocabulary, what's left genuinely unsupported by GlyLES itself (not just
unmapped) is documented above: `n_methyl`, `n_formyl`, `n_succinate`
(wrong chemistry via naive composition), `thio` (no sulfur bridge token
exists — a real grammar limitation), `pyrophosphate`/`triphosphate` (no
FG code exists for them at all), `pyruvate` (GlyLES only supports a
chemically-different single-position ester form, not the real
two-position cyclic acetal), `imino`/`dimethylamine`. A handful of rare
ones (`anhydro`, `lactone`, `diphospho_ethanolamine`, `hydroxymethyl`)
weren't attempted — low priority, and some would need dedicated
non-table code paths rather than a `SUBSTITUENT_MAP` entry (`anhydro` in
particular uses a completely different grammar rule, not bridge+FG).
**None of these are expected to come up often in real biological
glycans** — the 19 mapped entries cover every substituent type actually
seen in the test cases plus the other common ones (halogens, phosphate,
sulfate, common esters/amides). Treat this as "good enough to scale-test
against," not "100% complete" — the scale test (§2) will reveal whether
any of the unsupported ones actually matter in practice, which is more
useful than continuing to guess at rare cases with no real data to check
against.

Base-sugar identity (the bigger of the original two problems) is
effectively solved generically by `identity.identify()` — no comparable
per-sugar table is needed at all.

### 2. Scale testing — DONE this session, real number obtained, findings are mixed

Every step of this whole investigation has found that things which work
on clean/simple examples degrade on real, messy data (GlyLES's paper
claimed ~99%, real GlyTouCan data gave ~54%; the live production chain
looked fine until tested at scale and turned out to be 0%). This one is
no exception.

**Method** (apples-to-apples with the earlier ~54%/0% comparison in
`docs/SESSION_HANDOFF.md`): reused the *exact same* 100-row sample
(`random_state=7` over the user's GlycoSmos bulk export,
`~/Downloads/download.tsv`, same 100 GlyTouCan IDs as the earlier
GlyLES-direct-from-IUPAC vs. live-chain comparison). For each ID, fetched
its `wurcs` field from the same live `gtcid2seqs` API the pipeline
already calls (Context A's actual data source — not synthetic data), ran
it through `glypy.io.wurcs.loads()` → v2 `translate()` → `glyles.convert()`
(with a 15s `SIGALRM` timeout per the known hang risk), and categorized
every outcome.

**Result: 41/100 (41%) fully successful.**

|  | live chain (§ handoff doc) | GlyLES direct from bulk IUPAC (§ handoff doc) | WURCS→glypy→v2→GlyLES (this test) |
|---|---|---|---|
| success rate, same 100 IDs | 0% | 54% | **41%** |

This is a genuine improvement over the currently-broken live chain (0% →
41%) but **does not beat the existing GlyLES-direct-from-bulk-IUPAC
route** on this same sample, contrary to what the plan's "suggested order
of work" hoped for. See "Revised recommendation" below for what this
means for integration.

**Two more real bugs were found and fixed via this test** (both already
folded into the code above — the numbers here are *after* the fixes,
starting from an initial 6% before them):

1. **Anomer-sensitivity was the single biggest blocker (6% → 41% jump).**
   `identify()` failed to name a residue at all whenever its anomer was
   undefined (`x`) — which is *extremely* common in real WURCS,
   most fundamentally because **a released oligosaccharide's reducing end
   genuinely has no fixed anomeric configuration** (it mutarotates in
   solution) — this isn't a data-quality problem, it's real chemistry, and
   every multi-residue glycan hits it at least once. Fixed by retrying
   `identify()` with the candidate's anomer forced to each anomer actually
   used in glypy's own registry (beta covers the vast majority; alpha
   covers the `aMan`-style exception) until one matches — since a sugar's
   core identity (Glc vs. Gal vs. Man, ...) doesn't depend on its anomeric
   configuration, which is tracked separately anyway. See
   `_ANOMER_CANDIDATES`/`_base_sugar_name` in the code above.
2. **`UnknownPosition` (-1) was being emitted as literal text.** glypy
   represents an undefined linkage/substituent position with the int
   sentinel `-1`, and the first version of this code just `f"..."`ed it
   directly into the output string, producing garbage like `?1--1`
   instead of GlyLES's expected `?1-?`. Fixed with a small `_pos()` helper
   used everywhere a position is formatted. This fix was necessary but,
   on its own, didn't change the overall success count on this
   particular sample — see finding 4 below for why.

**Failure breakdown (59 non-successes, categorized)**:

- **31% "glyles\_empty"**: the translator produced well-formed,
  syntactically valid IUPAC-condensed text, but `glyles.convert()`
  returned an empty SMILES for it. Inspected these directly: **29 of 31
  have every single linkage position in the molecule marked unknown**
  (`?` for both parent and child position throughout, e.g.
  `Gal(?1-?)Glc2NAc(?1-?)Gal(?1-?)...`). This is not a translator defect
  — with literally no known attachment point anywhere in the graph,
  GlyLES cannot build any 3D structure, correctly. It's the same
  fundamental "real glycan data is often genuinely, irreducibly
  ambiguous" finding this whole investigation keeps running into (echoes
  the earlier finding that pervasive `?` wildcards cap GlyLES's
  direct-from-IUPAC success rate too) — just surfacing at the
  GlyLES-resolution step instead of the translation step this time.
- **26% "translate\_fail"**: 25 of 26 are `identify()` genuinely unable to
  name the base sugar at all, even after the anomer-retry fix — inspecting
  these, the residue's own composition string (e.g. `axxxxh`) has
  **undefined stereocenters baked into glypy's own parsed representation**
  (glypy logs `UserWarning: Cannot infer chirality from 'axxxxh-1x_1-5'`
  for exactly these), meaning the *source WURCS itself* records a
  "generic hexose" without committing to Glc/Gal/Man/etc. identity — again
  a real data-ambiguity ceiling, not a fixable bug. The remaining 1 is the
  known `n_methyl` substituent gap (documented above).
- **2% "translate\_exception"**: both are `glypy`'s own WURCS parser
  hitting notation it doesn't support (`ValueError` on an ambiguous
  intra-residue `4|6` substituent position, `KeyError` on an
  unrecognized substituent composition string `OC/2C`) — failures
  upstream of this translator entirely, in glypy's parser itself, not
  something this code can fix.

**Revised recommendation** (supersedes the plan's original framing that
this would clearly "beat" the direct-IUPAC baseline): the WURCS route is
a real, substantial win for **Context A specifically**, because Context
A's current live chain is confirmed broken at 0% — going to 41% offline
is a strict improvement with no downside. But it should be **layered
alongside, not instead of**, the existing bulk-IUPAC-file-based GlyLES
route where that data is available (54% on the same molecules) — e.g.
try direct-from-bulk-IUPAC first when the bulk file covers a given
GlyTouCan ID, fall back to WURCS→glypy→v2→GlyLES only when it doesn't
(recall from `docs/cognate_ligands_generation_plan.md`: the bulk file
only has 124,957 of 254,097 populated rows) or when the direct route
fails for that specific ID. Whether to also combine both routes'
*successes* when they disagree, or simply prefer whichever succeeds
first, is an open question for whoever does the integration work (§3) —
not investigated this session.

### 2b. Context B (PDB-native WURCS) scale test — DONE, much stronger result

The 41% number above was measured against a random sample of
**GlyTouCan-registry** WURCS (from the bulk export file, i.e. entries
depositors submitted directly to GlyTouCan, covering the huge diversity
of registered structures, many partially/ambiguously specified). Context
B's real input population is different: WURCS descriptors attached to
**bound carbohydrate entities in solved PDB structures**
(`_pdbx_entity_branch_descriptor`, `type="WURCS"`) — i.e. only structures
crystallographers actually resolved and deposited, which was hypothesized
to be more completely/unambiguously specified on average. This was worth
checking before recommending an approach for Context B, rather than
assuming the Context A number transfers.

**Method**: queried RCSB's search API for all PDB entries with at least
one branched entity (14,759 total, retrieved 10,000 via
`rcsb_entry_info.branched_entity_count > 0`), took a `random.seed(7)`
sample of 100, and pulled each entry's WURCS descriptor directly via
RCSB's GraphQL API (`https://data.rcsb.org/graphql`, querying
`branched_entities { pdbx_entity_branch_descriptor { descriptor type } }`
and taking the first `type: "WURCS"` entry) — this mirrors exactly what
`_pdbx_entity_branch_descriptor.descriptor` contains in the real mmCIF
files Context B's code reads, just fetched via API instead of parsing
mmCIF directly. All 100 sampled entries had a WURCS descriptor. Ran the
same glypy→v2→GlyLES chain (with the same 15s timeout guard) as the
Context A test.

**Result: 96/100 (96%)** — dramatically higher than the 41% seen on
GlyTouCan-registry data, confirming the hypothesis: PDB-deposited
carbohydrate structures are far more completely and unambiguously
specified than the general GlyTouCan registry population.

**One more real, generically-fixable bug found via this test** (already
folded into the code above): the first run of this test came back at
91%, with 5 of the 9 failures all being the same family — sucrose-type
disaccharides (`Glc(a1-2)Fructofuranose`) where `glyles.convert()`
returned empty SMILES. Root cause: glypy's named-residue registry has a
*separate* entry, `"Fructofuranose"`, distinct from `"Fru"`, for the
furanose ring form specifically; with `identify()`'s default
`ignore_ring=True`, both structurally match a real furanose-ring fructose
equally well, and `"Fructofuranose"` happened to win by registry list
order. GlyLES has no `"Fructofuranose"` token at all — it uses a generic
convention instead (append a bare `f` suffix to any base sugar name for
its furanose form, e.g. `Fruf`), the same convention GlyLES's own
`MonomerFactory.__getitem__` implements for *any* sugar capable of both
ring forms. Fixed generically, not as a Fru-specific patch: blacklisted
`"Fructofuranose"` so `identify()` falls through to the short `"Fru"`
form, then append `"f"` in `_base_sugar_name` whenever the *original*
node's `ring_type` is furanose. Verified via real GlyLES output
(`Glc(a1-2)Fruf` → correct sucrose SMILES) and by re-running the full
100-sample test, which went from 91% → 96% with zero regressions on the
other 91 successes.

**Remaining 4 failures, categorized** (all genuinely exotic/rare
molecules, not systemic issues):
- 2 are glypy's own WURCS parser explicitly rejecting notation it doesn't
  support: `WURCSFeatureNotSupported: Bridging MAPs are not supported`
  (a modification bridging two separate residues, e.g. a cyclic
  cross-link) and a `KeyError` on an unrecognized complex
  glycolipid-style substituent string (long branched acyl/lipid chain
  notation) — both upstream of this translator, in glypy's parser itself.
- 2 are `identify()` genuinely unable to name the base sugar: a
  fluorinated open-chain pentose (synthetic/labeled sugar analog, not a
  natural residue) and what's very likely a Δ4,5-unsaturated uronic acid
  (a known glycosaminoglycan-lyase-digestion product motif, e.g. from
  heparin/heparan sulfate or chondroitin processing) — neither is in
  glypy's or GlyLES's simple named-registry vocabulary; both are
  legitimately exotic, not something a general-purpose translator table
  should be expected to cover from one example.

**Conclusion for Context B**: unlike Context A, where the live chain is
confirmed broken (making the offline route's exact success rate almost
moot — anything beats 0%), Context B's live chain does actually work
today (confirmed on at least one real case earlier this investigation).
But at **96% offline**, fully removing the live-API dependency (and its
flagged reliability risks: unthrottled per-structure calls, no
retry/backoff, dependency on CSDB's HTML-scraped endpoint) is a clear,
low-risk win for Context B specifically — much more clear-cut than
Context A's mixed 41%-vs-54% picture. **Recommendation: make the offline
translator the primary path for Context B**, keeping the existing live
chain only as a fallback for the ~4% (and structurally-exotic edge cases
generally), rather than the reverse.

### 2c. Context A run against the REAL cognate-ligand master set (not a random sample) — DONE, best result yet

Both §2 and §2b tested random *populations* (GlyTouCan registry at large;
PDB-deposited structures at large) as a proxy for real usage. This test
instead derives and runs against the **actual master set** of KEGG glycan
compounds that would flow through Context A in production — i.e. glycan
codes that are genuine substrate/product participants in KEGG reactions
tied to an EC number (mirroring what
`kegg_reaction_enzyme_df_exploded.entities` computes in
`get_ec_information.py`), not just "any glycan somebody registered."

**Method** (approximates the pipeline's own derivation without running
its full per-EC live pipeline, which would take hours): pulled two bulk
KEGG link tables in two single requests, no per-ID loop
(`rest.kegg.jp/link/reaction/enzyme` → 9,711 unique EC-linked reactions;
`rest.kegg.jp/link/glycan/reaction` → 464 unique glycan-participant
reactions), intersected them to get **392 unique KEGG glycan codes
(`G#####`) that are genuine reaction participants of an EC-linked
reaction** — this is the real master-set glycan population, an order of
magnitude smaller and more specific than the raw 11,274-entry KEGG GLYCAN
catalog cited in `docs/cognate_ligands_generation_plan.md` Question 3
(that number was the *whole* catalog, not filtered to reaction
relevance). Caveat: this intersects against reactions with *any*
EC-enzyme link, not the exact `EC_substrate_codes`/`EC_product_codes`
fields the real code also folds in from enzyme records directly — close
enough to be representative, not a byte-exact reproduction of the
pipeline's own dataframe.

For each of the 392, fetched its real KEGG flat-file record
(`rest.kegg.jp/get/{id}`, throttled ~3/sec) and parsed its `DBLINKS`
section for a GlyTouCan cross-reference (this is exactly the
`extract_secondary_id`/`get_kegg_compound_record` step in the real code).
For each unique GlyTouCan ID found, called the same live `gtcid2seqs` API
the pipeline calls, then ran the full glypy→v2 translator→GlyLES chain.

**Results**:

| Stage | Count | % of 392 |
|---|---|---|
| Total master-set KEGG glycan codes | 392 | 100% |
| ...with a GlyTouCan cross-reference at all | 284 | 72.4% |
| ...of those, live API returned a `wurcs` value | 274 | (274/284 = 96.5%) |
| ...of those, live API returned a `glycoct` value | **0** | **0%** — re-confirms the Context A live chain is completely dead on real master-set data too, not just the earlier random sample |
| **Translator chain succeeds** | **253** | **64.5% of the full 392-code master set, or 89.1% of the 284 that have a GlyTouCan ID at all** |

**89.1% (of translatable candidates) is notably higher than both earlier
tests (41% on the general GlyTouCan registry, closer to but still below
96% on PDB-deposited structures)** — makes sense on reflection: KEGG's
own glycan compounds tend to be core, canonical, well-characterized
metabolic intermediates (N-glycan precursors, simple disaccharides,
sugar-nucleotide donors) rather than the long tail of complex/ambiguous
structures across the *entire* GlyTouCan registry that dragged the §2
number down.

**The 108 (27.6%) with no GlyTouCan cross-reference at all are a
separate, structural gap** — not a translator problem, and not fixable
by improving the translator, the bulk file, or the live chain, since none
of those approaches can work without a GlyTouCan ID to start from. This
is the same category of gap noted in `docs/cognate_ligands_generation_plan.md`'s
KEGG GLYCAN analysis (93.3% of the full catalog lacks external coverage)
— just now quantified for the master-set-relevant subset specifically.

**Checked whether layering in the bulk-IUPAC-file route (§2's hybrid
recommendation for Context A) would rescue any of the 26 translator
failures here**: merged the 274 attempted GlyTouCan IDs against the bulk
file — 249 of them are present in it. Of the 21 translator failures among
those 249, only 2 have bulk-file coverage, and **both of those 2 also
fail via the bulk-file→GlyLES route** (both are `MurNAc(...)(a1-` —
truncated/open-reducing-end fragments, the exact issue the earlier
GlyLES-direct investigation found needs explicit stripping before GlyLES
can parse it, which wasn't applied here). **So for this specific
population, the hybrid layering adds ~0 net benefit over the translator
alone** — a different conclusion than §2's general-population test, where
the two routes covered meaningfully different ground. Worth noting for
whoever does integration (§3): the earlier "always try bulk-file first"
recommendation is still sound as a general strategy, it just won't move
the needle much for this particular master-set-glycan population
specifically.

**One more real bug found, root-caused, and confirmed to be upstream of
this translator (in glypy itself, not fixable here)**: 21 of the 26
translator failures on this run are an identical
`TypeError: 'NoneType' object is not subscriptable`, always on WURCS
strings ending `~n` (e.g. `WURCS=2.0/1,1,1/[a2122h-1b_1-5]/1/a1-a3~n`) —
WURCS's **polymer/repeat-unit notation** ("this linkage repeats n
times"), used for generic homopolymer/backbone descriptions. This
notation is common in KEGG's glycan compounds specifically (many are
described as generic repeat structures, not one concrete molecule) but
didn't show up at all in the §2/§2b random-population tests. Traced the
full traceback: it originates inside **glypy's own**
`Glycan.__init__` → `canonicalize()` → `compare_residue_ordering()` →
`get_branch_from_link_label()`, which does `link.label[0]` on a link
whose `label` is `None` for this WURCS feature — a genuine glypy
limitation on repeat-unit WURCS, not a bug in this translator's code
(the exception happens before `translate()`'s own logic even runs). Left
uncaught-but-categorized (already correctly bucketed as
`translate_exception` in the scale-test harness) rather than worked
around, since fixing glypy's internals is out of scope here — worth
flagging to glypy upstream if this population's coverage matters enough
to chase further, but not a small fix.

### 3. Integration into the actual pipeline (not started)

Once the translator is solid and its real success rate is known, it needs
wiring into the actual codebase:

- **Context A** (`nextflow/bin/get_ec_information.py`, feeds
  `cognate_ligands_df.pkl`): replace or supplement the currently-broken
  `get_gtc_info`→CSDB chain with GlyTouCan-API-`wurcs` → glypy →
  translator → GlyLES.
- **Context B** (`nextflow/bin/process_all_pdb_contacts.py`, bound-entity
  glycans already observed in PDB structures — this path was tested and
  found to still work via its live chain, so this is an optional
  reliability upgrade here, not an urgent fix): could replace
  `get_glycoct_from_wurcs`→CSDB→CSDB with glypy → translator → GlyLES.
- Add `glypy` and `glyles` to `nextflow/envs/procoggraph.yaml`.
- Move the translator from scratchpad prototype into an actual module
  under `nextflow/bin/` (not done yet — the code above has only been run
  from a scratch venv, never placed in the repo or import-tested against
  the pipeline's actual environment/dependencies).
- Decide the fallback strategy: if the offline chain fails for a given
  structure, is there still value in falling back to the live API chain
  (Context B's, which was shown to still work for at least simple cases),
  or is offline-only acceptable given the live chain's own unresolved
  reliability issues (unthrottled per-item calls, no retry/backoff — see
  `docs/cognate_ligands_generation_plan.md`)?

## Suggested order of work

1. ~~Investigate GlyLES's/glypy's own internal vocabulary/object-model
   APIs before expanding `RESIDUE_MAP` by hand.~~ **Done this session** —
   led directly to v2.
2. ~~If switching to programmatic traversal of glypy's object model turns
   out to be straightforward, do that now, before scale-testing.~~ **Done
   this session** — v2 is built and validated on 4 cases.
3. ~~Finish substituent vocabulary coverage.~~ **Done this session** — 19
   entries, all individually verified against real GlyLES output; the
   remaining gaps are confirmed GlyLES-side limitations, not just
   unmapped cases (see §1 above).
4. ~~Run the scale test and get a real success-rate number.~~ **Done this
   session, for both contexts** — Context A: 41% on the same 100-ID
   sample used for the earlier 0%/54% comparison (§2). Context B: **96%**
   on a real sample of PDB-native WURCS descriptors (§2b) — much
   stronger, because PDB-deposited structures are far more completely
   specified than the general GlyTouCan registry population. Three real
   bugs found and fixed across both runs (see §2/§2b).
5. **Next concrete step**: decide the integration approach and wire it in
   (§3). The two contexts now have different clear recommendations:
   **Context B** — make the offline translator the *primary* path (96%,
   low-risk win, removes a live-API dependency that isn't broken but
   carries known reliability risk), live chain as fallback only.
   **Context A** — layer it with the existing bulk-file route: try
   direct-from-bulk-IUPAC first where the bulk file covers a GlyTouCan ID
   (54%), fall back to this translator otherwise (41%) — the live chain
   there is confirmed broken (0%) so isn't a fallback option at all.
   ~~Either way, this includes moving the code from prototype into an
   actual `nextflow/bin/` module.~~ **Done** —
   `nextflow/bin/wurcs_to_iupac.py`. Still not imported by either
   `get_ec_information.py` or `process_all_pdb_contacts.py`, and
   `glypy`/`glyles` aren't in `nextflow/envs/procoggraph.yaml` yet — the
   actual call-site wiring and the layering/fallback logic above are what
   remain.
6. Add a subprocess timeout wrapper around GlyLES calls given the
   hang-on-malformed-input behavior found during prototyping (the scale
   test used a 15s `SIGALRM` wrapper for this — reuse that approach in
   production code, adjusted for whatever concurrency model
   `nextflow/bin/get_ec_information.py` ends up using).

## Superseded: v1 (regex-over-text prototype, kept for history only)

The original prototype worked by regexing glypy's **text** IUPAC-extended
output, with a deliberately tiny 5-entry `RESIDUE_MAP`
(`Glcp2NAc→GlcNAc`, `Manp→Man`, `Galp→Gal`, `Glcp→Glc`, `Fucp→Fuc`). It
proved the same two cases (chitobiose, Man3GlcNAc2) work end-to-end and
that branch nesting is a token-level text transform, not a
graph-restructuring problem — both findings still hold and are re-proven
by v2. v2 fully replaces v1's implementation; v1's code is no longer
needed and is omitted here (see git history / earlier version of this
file if it's ever needed for reference).

## Where this fits in the bigger picture

This is a sub-thread of `docs/cognate_ligands_generation_plan.md` (see
its "Question 5" and "Question 6" sections for the full investigation
this grew out of, including the KEGG-side findings and the currently
broken production glycan path this is meant to fix). `docs/SESSION_HANDOFF.md`
has a running log of every finding from that investigation if more
context is needed than what's summarized here.
