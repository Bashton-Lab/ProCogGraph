#!/usr/bin/env python

"""
get_ec_information.py's cognate_ligands_df concat/groupby block is inline
script code (not a function), so - same approach as
test_get_ec_information_glycan_fallback.py - this replicates its exact
logic against small synthetic dataframes, to check the ligand_source
plumbing added for docs/cofactor_coverage_plan.md without needing the
full pipeline's upstream KEGG/Rhea machinery.

    python3 nextflow/bin/tests/test_get_ec_information_cofactor_merge.py
"""

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def make_frame(rows):
    df = pd.DataFrame(rows)
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol="smiles", molCol="ROMol")
    return df


class TestCofactorLigandSourceMerge(unittest.TestCase):

    def setUp(self):
        # one reaction-derived compound (unique), and one compound that
        # appears via BOTH a reaction source and the cofactor source (same
        # canonical SMILES, different EC/source) to check the union case.
        self.reaction_df = make_frame([
            {"entry": "1.1.1.1", "compound_id": "C00003", "compound_name": "NAD+",
             "smiles": "CC(=O)C", "ligand_db": "KEGG:C00003", "compound_reaction": "R00001"},
            {"entry": "2.2.2.2", "compound_id": "C00099", "compound_name": "some substrate",
             "smiles": "CCO", "ligand_db": "KEGG:C00099", "compound_reaction": "R00099"},
        ]).assign(ligand_source="reaction")

        self.cofactor_df = make_frame([
            # same compound as the first reaction row (CC(=O)C), but a
            # different EC, sourced from the cofactor path
            {"entry": "3.3.3.3", "compound_id": "NAD", "compound_name": "NAD",
             "smiles": "CC(=O)C", "ligand_db": "CofactorDB:4", "compound_reaction": ""},
        ]).assign(ligand_source="cofactor")

    def canon_smiles(self, x):
        try:
            return Chem.MolToSmiles(x, isomericSmiles=False)
        except Exception:
            return np.nan

    def run_merge_logic(self, reaction_df, cofactor_df):
        # exact logic from get_ec_information.py's cognate_ligands_df block
        cols = ["entry", "compound_id", "compound_name", "ROMol", "ligand_db", "compound_reaction", "ligand_source"]
        cognate_ligands_df = pd.concat([reaction_df[cols], cofactor_df[cols]])
        cognate_ligands_df = cognate_ligands_df.reset_index()

        cognate_ligands_df["canonical_smiles"] = cognate_ligands_df["ROMol"].map(
            lambda x: self.canon_smiles(x) if isinstance(x, Chem.rdchem.Mol) else np.nan
        )
        cognate_ligands_df_unique_smiles = cognate_ligands_df[
            ["canonical_smiles", "compound_name", "ligand_db", "compound_reaction", "ligand_source"]
        ].copy()
        cognate_ligands_df_unique_smiles["compound_reaction"] = cognate_ligands_df_unique_smiles["compound_reaction"].fillna("")
        cognate_ligands_df_unique_smiles = cognate_ligands_df_unique_smiles.groupby("canonical_smiles", dropna=False).agg(
            {"compound_name": set, "ligand_db": set, "compound_reaction": set, "ligand_source": set}
        ).reset_index()
        cognate_ligands_df_unique_smiles["ligand_db"] = cognate_ligands_df_unique_smiles.ligand_db.str.join("|")
        cognate_ligands_df_unique_smiles["compound_name"] = cognate_ligands_df_unique_smiles.compound_name.str.join("|")
        cognate_ligands_df_unique_smiles["compound_reaction"] = cognate_ligands_df_unique_smiles.compound_reaction.str.join("|").str.strip("|")
        cognate_ligands_df_unique_smiles["ligand_source"] = cognate_ligands_df_unique_smiles.ligand_source.apply(lambda x: "|".join(sorted(x)))
        cognate_ligands_df_unique_smiles = cognate_ligands_df_unique_smiles.reset_index(drop=True).reset_index()
        cognate_ligands_df_unique_smiles.rename(columns={"index": "uniqueID"}, inplace=True)

        cognate_ligands_df = cognate_ligands_df[["entry", "canonical_smiles"]].merge(
            cognate_ligands_df_unique_smiles, on="canonical_smiles", how="left"
        )
        return cognate_ligands_df.drop_duplicates()

    def test_reaction_only_compound_keeps_reaction_source(self):
        result = self.run_merge_logic(self.reaction_df, pd.DataFrame(columns=self.cofactor_df.columns))
        row = result.loc[result["entry"] == "2.2.2.2"].iloc[0]
        self.assertEqual(row["ligand_source"], "reaction")

    def test_compound_shared_between_reaction_and_cofactor_sources_gets_union_tag(self):
        result = self.run_merge_logic(self.reaction_df, self.cofactor_df)
        # both the reaction-path row (EC 1.1.1.1) and the cofactor-path row
        # (EC 3.3.3.3) point at the same canonical_smiles, so both entries
        # should carry the combined ligand_source tag.
        for ec in ["1.1.1.1", "3.3.3.3"]:
            row = result.loc[result["entry"] == ec].iloc[0]
            self.assertEqual(row["ligand_source"], "cofactor|reaction")

    def test_cofactor_only_entry_is_present_with_correct_ec(self):
        result = self.run_merge_logic(self.reaction_df, self.cofactor_df)
        self.assertIn("3.3.3.3", result["entry"].tolist())

    def test_existing_reaction_rows_are_unaffected_by_cofactor_addition(self):
        """docs/cofactor_coverage_plan.md's verification checklist: this
        must be a strict addition - the unrelated reaction-only compound's
        row should be identical whether or not any cofactor data is
        supplied."""
        without_cofactors = self.run_merge_logic(self.reaction_df, pd.DataFrame(columns=self.cofactor_df.columns))
        with_cofactors = self.run_merge_logic(self.reaction_df, self.cofactor_df)

        row_without = without_cofactors.loc[without_cofactors["entry"] == "2.2.2.2"].iloc[0]
        row_with = with_cofactors.loc[with_cofactors["entry"] == "2.2.2.2"].iloc[0]
        self.assertEqual(row_without["ligand_source"], row_with["ligand_source"])
        self.assertEqual(row_without["ligand_db"], row_with["ligand_db"])


if __name__ == "__main__":
    unittest.main()
