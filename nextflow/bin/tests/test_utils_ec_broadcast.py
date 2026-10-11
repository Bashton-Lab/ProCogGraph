#!/usr/bin/env python

"""
Tests for the EC-completeness classification and subsubclass-level
broadcast helpers added to utils.py for docs/cofactor_coverage_plan.md,
plus the get_chem_comp_descriptors CCD->SMILES resolver (moved here from
process_all_pdb_contacts.py, now reused by preprocess_cofactors.py).

    python3 nextflow/bin/tests/test_utils_ec_broadcast.py
"""

import sys
import unittest
from pathlib import Path

from gemmi import cif

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils import classify_ec_completeness, broadcast_subsubclass_ec, get_chem_comp_descriptors


class TestClassifyEcCompleteness(unittest.TestCase):

    def test_fully_resolved_ec_is_level_4(self):
        self.assertEqual(classify_ec_completeness("1.1.1.10"), 4)

    def test_subsubclass_level_wildcard_is_level_3(self):
        self.assertEqual(classify_ec_completeness("1.1.1.-"), 3)

    def test_subclass_level_wildcard_is_level_2(self):
        self.assertEqual(classify_ec_completeness("1.1.-.-"), 2)

    def test_class_level_wildcard_is_level_1(self):
        self.assertEqual(classify_ec_completeness("1.-.-.-"), 1)


class TestBroadcastSubsubclassEc(unittest.TestCase):

    def setUp(self):
        self.terminal_ec_list = ["1.1.1.1", "1.1.1.10", "1.1.2.1", "2.1.1.1"]

    def test_matches_only_same_subsubclass(self):
        matches = broadcast_subsubclass_ec("1.1.1.-", self.terminal_ec_list)
        self.assertCountEqual(matches, ["1.1.1.1", "1.1.1.10"])

    def test_no_matches_returns_empty_list(self):
        matches = broadcast_subsubclass_ec("3.3.3.-", self.terminal_ec_list)
        self.assertEqual(matches, [])

    def test_exact_ec_matches_itself_only(self):
        matches = broadcast_subsubclass_ec("1.1.1.1", self.terminal_ec_list)
        self.assertEqual(matches, ["1.1.1.1"])


class TestGetChemCompDescriptors(unittest.TestCase):
    """Reproduces the OpenEye-preference and invalid-SMILES-filtering
    branches with a small hand-written CCD-shaped CIF document, rather
    than needing the real (multi-GB) ccd.cif."""

    CCD_TEXT = """
data_XXX
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
XXX SMILES ACD 12.0 CCO
XXX SMILES 'OpenEye OEToolkits' 2.0.0 CCO
#
data_YYY
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
YYY SMILES ACD 12.0 not_a_valid_smiles(((
#
"""

    def setUp(self):
        self.ccd_doc = cif.Document()
        self.ccd_doc.parse_string(self.CCD_TEXT)

    def test_prefers_openeye_descriptor_when_present(self):
        result = get_chem_comp_descriptors(self.ccd_doc, ["XXX"])
        self.assertEqual(result["XXX"], "CCO")

    def test_invalid_smiles_resolves_to_none(self):
        result = get_chem_comp_descriptors(self.ccd_doc, ["YYY"])
        self.assertIsNone(result["YYY"])

    def test_absent_ligand_resolves_to_none(self):
        result = get_chem_comp_descriptors(self.ccd_doc, ["ZZZ"])
        self.assertIsNone(result["ZZZ"])


if __name__ == "__main__":
    unittest.main()
