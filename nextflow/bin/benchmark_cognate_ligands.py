#!/usr/bin/env python

"""
Compares two cognate_ligands_df.pkl files (e.g. a pre-cofactor-coverage
baseline vs. a new build) and reports the headline coverage numbers:
row/EC/ligand counts, which EC-ligand pairs are gained/lost, a rough
categorization of lost pairs, and the isCofactor tag distribution.

Written to check docs/cofactor_coverage_plan.md's real-world payoff
(see "Benchmark results" there for the run this script was built to
reproduce) - re-run it after any future change to the cognate-ligand
generation pipeline to see the actual before/after effect, not just
trust that the code changed.

Example usage:
    python3 benchmark_cognate_ligands.py \
        --old /path/to/old_cognate_ligands_df.pkl \
        --new /path/to/new_cognate_ligands_df.pkl \
        --report_out benchmark_report.txt
"""

import argparse
import re
import sys

import pandas as pd

WILDCARD_PATTERN = re.compile(r"\*")
MACROCYCLE_PATTERN = re.compile(r"\[N\+\]|\[Mg")


def load(path):
    df = pd.read_pickle(path)
    required = {"entry", "canonical_smiles", "compound_name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing expected column(s): {missing}")
    return df


def summarise_counts(df):
    return {
        "rows": len(df),
        "distinct_ecs": df["entry"].nunique(),
        "distinct_ligands": df["canonical_smiles"].nunique(),
    }


def pair_diff(old_df, new_df):
    old_pairs = set(zip(old_df["entry"], old_df["canonical_smiles"]))
    new_pairs = set(zip(new_df["entry"], new_df["canonical_smiles"]))
    return {
        "shared": old_pairs & new_pairs,
        "new_only": new_pairs - old_pairs,
        "old_only": old_pairs - new_pairs,
    }


def categorise_lost_pairs(old_df, old_only_pairs):
    """Rough heuristic split of pairs present in the old dataset but not
    the new one: wildcard-substituent partial structures, charge-separated
    porphyrin/macrocycle-like structures (both plausible RDKit-sanitization
    edge cases), vs. everything else (more likely ordinary upstream
    Rhea/KEGG reaction-equation revisions between the two builds - see the
    plan doc for a worked example: EC 1.1.1.96 listing NAD where it used
    to list NADH)."""
    lost_df = old_df[old_df.apply(lambda r: (r["entry"], r["canonical_smiles"]) in old_only_pairs, axis=1)]
    lost_structures = lost_df.drop_duplicates(subset="canonical_smiles")

    has_wildcard = lost_structures["canonical_smiles"].str.contains(WILDCARD_PATTERN, regex=True)
    has_macrocycle = lost_structures["canonical_smiles"].str.contains(MACROCYCLE_PATTERN, regex=True)

    return {
        "distinct_lost_structures": len(lost_structures),
        "wildcard_substituent": int(has_wildcard.sum()),
        "charge_separated_macrocycle": int((has_macrocycle & ~has_wildcard).sum()),
        "other": int((~has_wildcard & ~has_macrocycle).sum()),
    }


def format_report(old_path, new_path, old_counts, new_counts, diff, lost_ecs, lost_categories, old_cofactor, new_cofactor):
    lines = []
    lines.append(f"Cognate ligand dataset benchmark")
    lines.append(f"  old: {old_path}")
    lines.append(f"  new: {new_path}")
    lines.append("")
    lines.append(f"{'metric':<28}{'old':>12}{'new':>12}{'change':>12}")
    for key, label in [("rows", "EC-ligand pairs"), ("distinct_ecs", "distinct ECs"), ("distinct_ligands", "distinct ligands")]:
        old_v, new_v = old_counts[key], new_counts[key]
        change = f"{new_v - old_v:+d}"
        lines.append(f"{label:<28}{old_v:>12}{new_v:>12}{change:>12}")
    lines.append("")
    lines.append(f"shared EC-ligand pairs: {len(diff['shared'])}")
    lines.append(f"new-only pairs (gained): {len(diff['new_only'])}")
    lines.append(f"old-only pairs (lost): {len(diff['old_only'])}")
    lines.append("")
    lines.append(f"ECs lost entirely (in old, absent from new): {len(lost_ecs)}")
    if lost_ecs:
        lines.append(f"  {', '.join(lost_ecs)}")
    lines.append("")
    lines.append("Lost-pair categorisation (distinct structures, heuristic):")
    lines.append(f"  distinct lost structures: {lost_categories['distinct_lost_structures']}")
    lines.append(f"  wildcard-substituent (partial structures): {lost_categories['wildcard_substituent']}")
    lines.append(f"  charge-separated macrocycle (e.g. porphyrins): {lost_categories['charge_separated_macrocycle']}")
    lines.append(f"  other (likely upstream Rhea/KEGG data drift): {lost_categories['other']}")
    lines.append("")
    if old_cofactor is not None and new_cofactor is not None:
        lines.append("isCofactor tag distribution:")
        lines.append(f"{'':<28}{'old':>12}{'new':>12}")
        all_tags = sorted(set(old_cofactor.index) | set(new_cofactor.index))
        for tag in all_tags:
            lines.append(f"{tag:<28}{old_cofactor.get(tag, 0):>12}{new_cofactor.get(tag, 0):>12}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Benchmark two cognate_ligands_df.pkl files against each other")
    parser.add_argument("--old", required=True, help="Path to the baseline cognate_ligands_df.pkl")
    parser.add_argument("--new", required=True, help="Path to the new cognate_ligands_df.pkl")
    parser.add_argument("--report_out", default=None, help="Optional path to also write the report to a file")
    args = parser.parse_args()

    old_df = load(args.old)
    new_df = load(args.new)

    old_counts = summarise_counts(old_df)
    new_counts = summarise_counts(new_df)
    diff = pair_diff(old_df, new_df)

    lost_ecs = sorted(set(old_df["entry"]) - set(new_df["entry"]))
    lost_categories = categorise_lost_pairs(old_df, diff["old_only"])
    lost_categories["lost_ecs"] = lost_ecs

    old_cofactor = old_df["isCofactor"].value_counts() if "isCofactor" in old_df.columns else None
    new_cofactor = new_df["isCofactor"].value_counts() if "isCofactor" in new_df.columns else None

    report = format_report(args.old, args.new, old_counts, new_counts, diff, lost_ecs, lost_categories, old_cofactor, new_cofactor)
    print(report)

    if args.report_out:
        with open(args.report_out, "w") as fh:
            fh.write(report + "\n")
        print(f"\nReport written to {args.report_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
