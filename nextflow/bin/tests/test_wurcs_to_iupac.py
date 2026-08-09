#!/usr/bin/env python

"""
Small, offline-only regression tests for wurcs_to_iupac.translate().
Locks in the cases validated during development (see
docs/iupac_translator_plan.md) so future changes can't silently regress
them. No network access required - run with:

    python3 nextflow/bin/tests/test_wurcs_to_iupac.py
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from wurcs_to_iupac import translate, TranslationError, SUBSTITUENT_MAP, _substituent_token
from glypy.structure.substituent import Substituent


class TestTranslate(unittest.TestCase):

    def test_simple_disaccharide(self):
        # chitobiose core, pulled from real PDB entry 6VXX
        wurcs = "WURCS=2.0/1,2,1/[a2122h-1b_1-5_2*NCC/3=O]/1-1/a4-b1"
        self.assertEqual(translate(wurcs), "Glc2NAc(b1-4)Glc2NAc")

    def test_branched_n_glycan_core(self):
        # Man3GlcNAc2 - the standard N-glycosylation core
        wurcs = (
            "WURCS=2.0/3,5,4/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1a_1-5][a1122h-1a_1-5]"
            "/1-1-2-3-3/a4-b1_b4-c1_c3-d1_c6-e1"
        )
        self.assertEqual(
            translate(wurcs),
            "Man(a1-6)[Man(a1-3)]Man(a1-4)Glc2NAc(b1-4)Glc2NAc",
        )

    def test_furanose_ring_suffix(self):
        # sucrose - regression test for the "Fructofuranose" vs "Fruf"
        # naming quirk found during PDB-scale testing
        wurcs = "WURCS=2.0/2,2,1/[ha122h-2b_2-5][a2122h-1a_1-5]/1-2/a2-b1"
        self.assertEqual(translate(wurcs), "Glc(a1-2)Fruf")

    def test_undefined_anomer_reducing_end(self):
        # a free reducing end has no fixed anomeric configuration in real
        # WURCS data - regression test for the anomer-override retry in
        # _base_sugar_name (this used to raise TranslationError entirely)
        wurcs = "WURCS=2.0/1,1,1/[a2122h-1x_1-5]/1/"
        result = translate(wurcs)
        self.assertTrue(result.startswith("Glc"))

    def test_unsupported_substituent_raises(self):
        # n_methyl is confirmed-unsupported (composes to the wrong
        # chemistry in GlyLES, see docs/iupac_translator_plan.md) and
        # must not be in the map, and must fail loudly rather than guess
        self.assertNotIn("n_methyl", SUBSTITUENT_MAP)
        with self.assertRaises(TranslationError):
            _substituent_token(2, Substituent("n_methyl"))

    def test_malformed_wurcs_does_not_silently_succeed(self):
        with self.assertRaises(Exception):
            translate("not a wurcs string")


if __name__ == "__main__":
    unittest.main()
