#!/usr/bin/env python

"""
Translate a WURCS glycan descriptor to GlyLES-compatible IUPAC-condensed
notation, via glypy's structured WURCS parse (not glypy's own text IUPAC
export, which is a different, incompatible dialect).

Background and full validation history: docs/iupac_translator_plan.md.

Benchmark comparison (glypy -> translate() -> glyles.convert(), vs. the
existing live-API chains, not re-run as part of this commit - see the plan
doc for full methodology):

Context A (cognate ligand master set, get_ec_information.py, WURCS sourced
from GlyTouCan's gtcid2seqs API):
  - existing live chain (glycoct via gtcid2seqs -> CSDB):        0%
    (GlyTouCan's API stopped returning `glycoct`; confirmed dead
    on every sample tested, including the real 392-glycan master set)
  - bulk GlycoSmos IUPAC-condensed export -> glyles directly:   54%
    (100-sample random GlyTouCan registry ID)
  - this module (glypy -> translate() -> glyles):               41%
    (same 100-sample random GlyTouCan registry ID)
  - this module, against the REAL 392-glycan cognate-ligand
    master set specifically (not a random sample):             89.1%
    (of the 284/392 that have a GlyTouCan ID at all; 64.5% of
    the full 392 including the 108 with no GlyTouCan ID, which
    no approach here can resolve)
  Recommendation: layer bulk-file-first, this module as fallback,
  for the general population - but for the real master set specifically,
  this module alone already covers about as much as the hybrid does.

Context B (bound-entity glycans in solved PDB structures,
process_all_pdb_contacts.py, WURCS sourced directly from
_pdbx_entity_branch_descriptor):
  - existing live chain (GlycoSmos wurcs2glycoct -> CSDB -> CSDB): works
    today, but unthrottled/no-retry, and depends on an HTML-scraped
    CSDB endpoint
  - this module:                                                  96%
    (100-sample of real PDB entries with branched entities, pulled
    live via RCSB's search + GraphQL APIs)
  Recommendation: this module as the primary path, live chain as
  fallback only.

Neither integration (replacing the relevant calls in get_ec_information.py
/ process_all_pdb_contacts.py) has been wired in yet - this module is not
yet imported anywhere in the pipeline.
"""

import glypy
from glypy.io import wurcs as wurcs_io
from glypy.io.nomenclature import identity
from glypy.structure.constants import UnknownPosition


class TranslationError(Exception):
    pass


def _pos(value):
    """glypy's int sentinel for an undefined position (-1) has to become
    GlyLES's '?' token, not its own literal text - emitting -1 directly
    produces unparseable linkage notation like '?1--1'."""
    return "?" if value == UnknownPosition else str(value)


# "Fructofuranose" is a separate glypy registry entry from "Fru" for the
# furanose ring form specifically. identify()'s default ignore_ring=True
# makes both match a real furanose-ring fructose equally well, and
# "Fructofuranose" wins by registry list order - but GlyLES has no such
# token, it appends a bare "f" suffix to the base name instead (handled
# generically below via node.ring_type, not just for fructose).
IDENTIFY_BLACKLIST = {"Pen", "Hex", "Hep", "Oct", "Non", "Fructofuranose"}

# glypy's registry has exactly one entry ("aMan") whose name bakes in a
# specific anomer rather than being anomer-generic. Anomer is tracked
# separately via node.anomer at the linkage-emission step regardless, so
# strip it back off here.
ANOMER_PREFIXED_ALIASES = {"aMan": "Man", "bMan": "Man"}

# glypy substituent .name -> (bridge_letter_or_None, FG_code_or_None).
# Every entry verified against real GlyLES output, not derived from the
# grammar shape - GlyLES's FG codes are hand-written SMILES fragments that
# already include their own attachment atom, so bridge+FG composition is
# only correct for some combinations (e.g. "n_methyl"/"n_formyl" compose
# to the wrong chemistry and are deliberately left out - see the plan doc
# for the rest of glypy's ~35-entry substituent vocabulary and which of
# the remainder are confirmed unsupported by GlyLES itself vs. just not
# attempted).
SUBSTITUENT_MAP = {
    "n_acetyl": ("N", "Ac"),
    "n_glycolyl": ("N", "Gc"),
    "n_sulfate": ("N", "S"),
    "n_amidino": ("N", "Am"),
    "sulfate": (None, "S"),
    "phosphate": (None, "P"),
    "methyl": (None, "Me"),
    "acetyl": (None, "Ac"),
    "amino": ("N", None),
    "glycolyl": (None, "Gc"),
    "formyl": (None, "Fo"),
    "ethyl": (None, "Et"),
    "ethanolamine": (None, "Etn"),
    "phospho_ethanolamine": ("P", "Etn"),
    "phospho_choline": ("P", "Cho"),
    "succinate": (None, "Suc"),
    "bromo": (None, "Br"),
    "chloro": (None, "Cl"),
    "fluoro": (None, "F"),
    "iodo": (None, "I"),
}

# identify() is anomer-sensitive (matches a residue's anomer against each
# registry candidate's own fixed anomer). Real WURCS very commonly has
# undefined anomer - most fundamentally because a released
# oligosaccharide's reducing end genuinely has no fixed anomeric
# configuration (mutarotates in solution), not a data-quality issue. Base
# sugar identity doesn't depend on anomeric configuration (tracked
# separately via node.anomer regardless), so retry with anomer forced to
# each anomer glypy's registry actually uses until one matches.
_ANOMER_CANDIDATES = [None, glypy.monosaccharides["Glc"].anomer, glypy.monosaccharides["Fuc"].anomer]


def _base_sugar_name(node):
    bare = node.clone()
    for pos, sub in list(bare.substituents()):
        bare.drop_substituent(pos, sub)
    last_error = None
    for anomer_override in _ANOMER_CANDIDATES:
        candidate = bare.clone()
        if anomer_override is not None:
            candidate.anomer = anomer_override
        try:
            name = identity.identify(candidate, blacklist=IDENTIFY_BLACKLIST)
            name = ANOMER_PREFIXED_ALIASES.get(name, name)
            # GlyLES marks a furanose ring form with a trailing "f" on an
            # otherwise pyranose-default base name, rather than having
            # separate named tokens per ring form.
            if node.ring_type is not None and node.ring_type.name == "furanose":
                name += "f"
            return name
        except identity.IdentifyException as e:
            last_error = e
    raise TranslationError(f"unrecognized base sugar for residue {node.id}: {last_error}") from last_error


def _substituent_token(position, substituent):
    name = substituent.name
    if name not in SUBSTITUENT_MAP:
        raise TranslationError(f"no GlyLES mapping for substituent '{name}' (position {position})")
    bridge, fg = SUBSTITUENT_MAP[name]
    # a bare bridge atom (fg=None, e.g. "amino") is itself a valid fgi in
    # glyles' grammar (`fgi: COUNT | bridge | FG`) - verified via Glc2N.
    return f"{_pos(position)}{bridge or ''}{fg or ''}"


def _residue_token(node):
    base = _base_sugar_name(node)
    subs = "".join(
        _substituent_token(pos, sub)
        for pos, sub in sorted(node.substituents(), key=lambda ps: str(ps[0]))
    )
    return base + subs


def _emit(node):
    """Recursively emit GlyLES IUPAC-condensed text for the subtree rooted
    at `node`, non-reducing-end-first (matching GlyLES's own bracket
    convention)."""
    # node.links holds both the link back to this node's own parent (if
    # any) and links to its children, keyed by position on this node -
    # only descend where this node is the parent side. Sort by position as
    # a string: positions can be ints or the unknown-position sentinel,
    # and real ambiguous WURCS data can give multiple children the same
    # (undefined) parent position, which made a naive sort by (position,
    # Link) tuples crash comparing two Link objects on ties.
    children = [
        (parent_pos, link) for parent_pos, link in sorted(node.links.items(), key=lambda kv: str(kv[0]))
        if link.parent.id == node.id
    ]
    branch_parts = []
    main_child = None
    # GlyLES convention: the last-positioned child continues the main
    # chain unbracketed; all others are bracketed side branches.
    for parent_pos, link in children:
        main_child = (link.child, link)
    for parent_pos, link in children:
        if main_child is not None and link.child.id == main_child[0].id:
            continue
        branch_parts.append(f"[{_emit_edge(link.child, link)}]")

    own = _residue_token(node)
    suffix = "".join(branch_parts)
    if main_child is not None:
        child, link = main_child
        return f"{_emit_edge(child, link)}{suffix}{own}"
    return f"{suffix}{own}"


def _emit_edge(child, link):
    subtree = _emit(child)
    anomer = "a" if child.anomer.name.startswith("alpha") else (
        "b" if child.anomer.name.startswith("beta") else "?"
    )
    return f"{subtree}({anomer}{_pos(link.child_position)}-{_pos(link.parent_position)})"


def translate(wurcs_string):
    """
    Translate a WURCS glycan descriptor to GlyLES-compatible IUPAC-condensed
    notation.

    Args:
        wurcs_string (str): a WURCS 2.0 descriptor, e.g. as returned by
            GlyTouCan's gtcid2seqs API or a PDB entry's
            _pdbx_entity_branch_descriptor.

    Returns:
        str: IUPAC-condensed notation suitable for glyles.convert().

    Raises:
        TranslationError: a residue's base sugar or a substituent has no
            known GlyLES-compatible encoding.
    """
    glycan = wurcs_io.loads(wurcs_string)
    return _emit(glycan.root)
