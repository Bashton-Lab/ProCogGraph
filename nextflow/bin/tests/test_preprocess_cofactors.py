#!/usr/bin/env python

"""
Tests for preprocess_cofactors.py, the new cofactor-EC association
preprocessing script for docs/cofactor_coverage_plan.md. No network
access or real reference data required - all inputs are small synthetic
fixtures built in-memory.

    python3 nextflow/bin/tests/test_preprocess_cofactors.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from gemmi import cif

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from preprocess_cofactors import (
    load_cofactor_db_table,
    parse_uniprot_cofactor_df,
    classify_and_split_ec_rows,
    resolve_chebi_smiles_map,
    build_cofactor_ligands_df,
)


class TestParseUniprotCofactorDf(unittest.TestCase):

    def test_structured_row_with_single_ec_and_chebi(self):
        raw = pd.DataFrame([
            {"Entry": "P00001", "EC number": "1.1.1.1",
             "Cofactor": "COFACTOR: Name=NAD(+); Xref=ChEBI:CHEBI:57540;"},
        ])
        result = parse_uniprot_cofactor_df(raw)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["ec"], "1.1.1.1")
        self.assertEqual(result.iloc[0]["chebi_id"], 57540)

    def test_free_text_only_row_is_dropped(self):
        # e.g. chlorophyll's real photosystem entries - Note= with no Xref=
        raw = pd.DataFrame([
            {"Entry": "P56766", "EC number": "1.97.1.12",
             "Cofactor": "COFACTOR: Note=P700 is a chlorophyll a/chlorophyll a' dimer."},
        ])
        result = parse_uniprot_cofactor_df(raw)
        self.assertTrue(result.empty)

    def test_multiple_ecs_on_one_row_are_exploded(self):
        raw = pd.DataFrame([
            {"Entry": "P00002", "EC number": "1.1.1.1; 1.1.1.2",
             "Cofactor": "COFACTOR: Name=Zn(2+); Xref=ChEBI:CHEBI:29105;"},
        ])
        result = parse_uniprot_cofactor_df(raw)
        self.assertCountEqual(result["ec"].tolist(), ["1.1.1.1", "1.1.1.2"])

    def test_multiple_chebi_xrefs_on_one_row_are_exploded(self):
        raw = pd.DataFrame([
            {"Entry": "P00003", "EC number": "1.1.1.1",
             "Cofactor": "COFACTOR: Name=[4Fe-4S]; Xref=ChEBI:CHEBI:49883; "
                         "Name=[2Fe-2S]; Xref=ChEBI:CHEBI:49601;"},
        ])
        result = parse_uniprot_cofactor_df(raw)
        self.assertCountEqual(result["chebi_id"].tolist(), [49883, 49601])

    def test_row_with_no_ec_is_dropped(self):
        raw = pd.DataFrame([
            {"Entry": "P00004", "EC number": float("nan"),
             "Cofactor": "COFACTOR: Name=NAD(+); Xref=ChEBI:CHEBI:57540;"},
        ])
        result = parse_uniprot_cofactor_df(raw)
        self.assertTrue(result.empty)

    def test_duplicate_ec_chebi_pairs_are_deduplicated(self):
        raw = pd.DataFrame([
            {"Entry": "P00005", "EC number": "1.1.1.1",
             "Cofactor": "COFACTOR: Name=NAD(+); Xref=ChEBI:CHEBI:57540;"},
            {"Entry": "P00006", "EC number": "1.1.1.1",
             "Cofactor": "COFACTOR: Name=NAD(+); Xref=ChEBI:CHEBI:57540;"},
        ])
        result = parse_uniprot_cofactor_df(raw)
        self.assertEqual(len(result), 1)


class TestClassifyAndSplitEcRows(unittest.TestCase):

    def setUp(self):
        self.terminal_ec_list = ["1.1.1.1", "1.1.1.10", "1.1.2.1", "2.1.1.1"]

    def test_exact_ec_passes_through_unchanged(self):
        df = pd.DataFrame([{"ec": "1.1.1.1", "chebi_id": 1}])
        kept, dropped = classify_and_split_ec_rows(df, self.terminal_ec_list)
        self.assertEqual(kept["entry"].tolist(), ["1.1.1.1"])
        self.assertTrue(dropped.empty)

    def test_subsubclass_wildcard_broadcasts_to_matching_terminal_ecs_only(self):
        df = pd.DataFrame([{"ec": "1.1.1.-", "chebi_id": 1}])
        kept, dropped = classify_and_split_ec_rows(df, self.terminal_ec_list)
        self.assertCountEqual(kept["entry"].tolist(), ["1.1.1.1", "1.1.1.10"])
        self.assertTrue(dropped.empty)

    def test_class_level_wildcard_is_dropped_not_broadcast(self):
        """Regression test for the bug caught during planning: naively
        broadcasting a class-level wildcard (e.g. "1.-.-.-") would match
        literally every terminal EC in that class, incorrectly claiming
        structurally unrelated enzymes share a cofactor. This must never
        reach `kept` - it belongs in `dropped` instead."""
        df = pd.DataFrame([{"ec": "1.-.-.-", "chebi_id": 1}])
        kept, dropped = classify_and_split_ec_rows(df, self.terminal_ec_list)
        self.assertTrue(kept.empty)
        self.assertEqual(len(dropped), 1)
        self.assertIn("too coarse", dropped.iloc[0]["reason"])

    def test_subclass_level_wildcard_is_dropped_not_broadcast(self):
        df = pd.DataFrame([{"ec": "1.1.-.-", "chebi_id": 1}])
        kept, dropped = classify_and_split_ec_rows(df, self.terminal_ec_list)
        self.assertTrue(kept.empty)
        self.assertEqual(len(dropped), 1)

    def test_subsubclass_wildcard_with_no_match_is_logged_as_dropped(self):
        df = pd.DataFrame([{"ec": "9.9.9.-", "chebi_id": 1}])
        kept, dropped = classify_and_split_ec_rows(df, self.terminal_ec_list)
        self.assertTrue(kept.empty)
        self.assertEqual(len(dropped), 1)

    def test_mixed_rows_only_safe_ones_survive(self):
        df = pd.DataFrame([
            {"ec": "1.1.1.1", "chebi_id": 1},   # exact -> kept
            {"ec": "1.1.1.-", "chebi_id": 2},   # subsubclass -> broadcast
            {"ec": "1.-.-.-", "chebi_id": 3},   # class-level -> dropped
        ])
        kept, dropped = classify_and_split_ec_rows(df, self.terminal_ec_list)
        # exact row (1 entry) + broadcast row (2 entries: 1.1.1.1, 1.1.1.10)
        self.assertEqual(len(kept), 3)
        self.assertEqual(len(dropped), 1)
        self.assertEqual(dropped.iloc[0]["chebi_id"], 3)


class TestLoadCofactorDbTable(unittest.TestCase):

    def test_merges_ec_csv_with_representative_ccd_from_details_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = Path(tmpdir) / "cofactor_ec.csv"
            csv_path.write_text('"EC_NO","COFACTOR_ID","SOURCE"\n"1.1.1.1",4,cofactor_db_2010\n"1.1.1.2",4,brenda\n')

            json_path = Path(tmpdir) / "cofactors_details.json"
            json_path.write_text('[{"template": "NAD2", "id": 4, "representative": "NAD", "threshold": 0.68}]')

            result = load_cofactor_db_table(str(csv_path), str(json_path))
            self.assertEqual(len(result), 2)
            self.assertTrue((result["representative_ccd"] == "NAD").all())
            self.assertCountEqual(result["ec"].tolist(), ["1.1.1.1", "1.1.1.2"])

    def test_cofactor_class_missing_from_details_json_is_dropped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = Path(tmpdir) / "cofactor_ec.csv"
            csv_path.write_text('"EC_NO","COFACTOR_ID","SOURCE"\n"1.1.1.1",99,cofactor_db_2010\n')

            json_path = Path(tmpdir) / "cofactors_details.json"
            json_path.write_text('[{"template": "NAD2", "id": 4, "representative": "NAD", "threshold": 0.68}]')

            result = load_cofactor_db_table(str(csv_path), str(json_path))
            self.assertTrue(result.empty)


class TestResolveChebiSmilesMap(unittest.TestCase):

    def test_resolves_known_chebi_id(self):
        structures_df = pd.DataFrame({"compound_id": [57540, 29105], "smiles": ["C1=CC...", "[Zn+2]"]})
        names_df = pd.DataFrame({"compound_id": [57540, 29105], "name": ["NAD(+)", "zinc(2+)"]})
        result = resolve_chebi_smiles_map([57540], structures_df, names_df)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["smiles"], "C1=CC...")
        self.assertEqual(result.iloc[0]["name"], "NAD(+)")

    def test_unknown_chebi_id_is_silently_excluded(self):
        structures_df = pd.DataFrame({"compound_id": [57540], "smiles": ["C1=CC..."]})
        names_df = pd.DataFrame({"compound_id": [57540], "name": ["NAD(+)"]})
        result = resolve_chebi_smiles_map([99999999], structures_df, names_df)
        self.assertTrue(result.empty)


class TestBuildCofactorLigandsDf(unittest.TestCase):
    """End-to-end integration across both sources with small synthetic
    fixtures - no real reference data or network access."""

    CCD_TEXT = """
data_NAD
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
NAD SMILES 'OpenEye OEToolkits' 2.0.0 CC(=O)C
#
"""

    def setUp(self):
        self.ccd_doc = cif.Document()
        self.ccd_doc.parse_string(self.CCD_TEXT)
        self.terminal_ec_list = ["1.1.1.1", "1.1.1.10", "2.2.2.2"]

        self.cofactor_db_df = pd.DataFrame([
            # subsubclass wildcard -> should broadcast to 1.1.1.1 and 1.1.1.10
            {"ec": "1.1.1.-", "cofactor_id": 4, "representative_ccd": "NAD", "source": "cofactor_db_2010"},
        ])
        self.uniprot_df = pd.DataFrame([
            {"ec": "2.2.2.2", "chebi_id": 29105},
        ])
        self.chebi_structures_df = pd.DataFrame({"compound_id": [29105], "smiles": ["[Zn+2]"]})
        self.chebi_names_df = pd.DataFrame({"compound_id": [29105], "name": ["zinc(2+)"]})

    def test_output_shape_and_columns(self):
        result, dropped = build_cofactor_ligands_df(
            self.cofactor_db_df, self.uniprot_df, self.terminal_ec_list,
            self.ccd_doc, self.chebi_structures_df, self.chebi_names_df,
        )
        expected_cols = {"entry", "compound_id", "compound_name", "ROMol", "ligand_db", "compound_reaction", "ligand_source"}
        self.assertEqual(set(result.columns), expected_cols)
        self.assertTrue((result["ligand_source"] == "cofactor").all())
        self.assertTrue((result["compound_reaction"] == "").all())

    def test_cofactor_db_row_broadcasts_to_both_terminal_ecs(self):
        result, _ = build_cofactor_ligands_df(
            self.cofactor_db_df, self.uniprot_df, self.terminal_ec_list,
            self.ccd_doc, self.chebi_structures_df, self.chebi_names_df,
        )
        nad_entries = result.loc[result["compound_id"] == "NAD", "entry"].tolist()
        self.assertCountEqual(nad_entries, ["1.1.1.1", "1.1.1.10"])

    def test_uniprot_row_resolves_via_chebi(self):
        result, _ = build_cofactor_ligands_df(
            self.cofactor_db_df, self.uniprot_df, self.terminal_ec_list,
            self.ccd_doc, self.chebi_structures_df, self.chebi_names_df,
        )
        zinc_rows = result.loc[result["compound_id"] == "CHEBI:29105"]
        self.assertEqual(len(zinc_rows), 1)
        self.assertEqual(zinc_rows.iloc[0]["entry"], "2.2.2.2")
        self.assertEqual(zinc_rows.iloc[0]["compound_name"], "zinc(2+)")

    def test_no_wildcard_ec_ever_appears_in_entry_column(self):
        """Hard invariant the downstream get_pdb_parity.py join depends on
        (docs/cofactor_coverage_plan.md's verification checklist)."""
        result, _ = build_cofactor_ligands_df(
            self.cofactor_db_df, self.uniprot_df, self.terminal_ec_list,
            self.ccd_doc, self.chebi_structures_df, self.chebi_names_df,
        )
        self.assertFalse(result["entry"].str.contains("-").any())

    def test_unresolvable_ccd_code_is_excluded_not_erroring(self):
        cofactor_db_df = pd.DataFrame([
            {"ec": "1.1.1.1", "cofactor_id": 999, "representative_ccd": "NOPE", "source": "cofactor_db_2010"},
        ])
        result, _ = build_cofactor_ligands_df(
            cofactor_db_df, pd.DataFrame(columns=["ec", "chebi_id"]), self.terminal_ec_list,
            self.ccd_doc, self.chebi_structures_df, self.chebi_names_df,
        )
        self.assertTrue(result.empty)


if __name__ == "__main__":
    unittest.main()
