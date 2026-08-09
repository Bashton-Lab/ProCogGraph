#!/usr/bin/env python

"""
Small tests for the offline glycan-to-SMILES helper added to utils.py,
and for how process_all_pdb_contacts.get_sugar_smiles_from_wurcs layers
it in front of the existing live GlycoSmos/CSDB chain. No network access
required (the live-chain calls are monkeypatched, not actually made).

    python3 nextflow/bin/tests/test_utils_glycan_offline.py
"""

import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils import get_smiles_from_wurcs_offline
import process_all_pdb_contacts as pac

CHITOBIOSE_WURCS = "WURCS=2.0/1,2,1/[a2122h-1b_1-5_2*NCC/3=O]/1-1/a4-b1"


class TestGetSmilesFromWurcsOffline(unittest.TestCase):

    def test_valid_wurcs_returns_smiles(self):
        smiles = get_smiles_from_wurcs_offline(CHITOBIOSE_WURCS)
        self.assertIsInstance(smiles, str)
        self.assertNotEqual(smiles, "")

    def test_none_returns_nan(self):
        self.assertTrue(pd.isna(get_smiles_from_wurcs_offline(None)))

    def test_nan_returns_nan(self):
        self.assertTrue(pd.isna(get_smiles_from_wurcs_offline(np.nan)))

    def test_malformed_wurcs_returns_nan_not_exception(self):
        # a production pipeline can't have a bad WURCS string crash the
        # whole run - failures must degrade to nan, same contract as the
        # existing get_glycoct_from_wurcs/get_csdb_from_glycoct/
        # get_smiles_from_csdb functions in utils.py
        self.assertTrue(pd.isna(get_smiles_from_wurcs_offline("not a wurcs string")))


class TestContextBFallbackWiring(unittest.TestCase):
    """
    process_all_pdb_contacts.get_sugar_smiles_from_wurcs should try the
    offline route first and only fall back to the live chain
    (get_glycoct_from_wurcs -> get_csdb_from_glycoct -> get_smiles_from_csdb)
    when the offline route fails - per the Context B recommendation in
    docs/iupac_translator_plan.md (offline primary, live chain fallback).
    """

    def _empty_cache(self, columns):
        return pd.DataFrame(columns=columns)

    def test_offline_success_skips_live_chain_entirely(self):
        with mock.patch.object(pac, "get_glycoct_from_wurcs") as mock_glycoct:
            sugar_smiles, *_ = pac.get_sugar_smiles_from_wurcs(
                [CHITOBIOSE_WURCS],
                self._empty_cache(["glycoct", "csdb"]),
                self._empty_cache(["csdb", "descriptor"]),
                self._empty_cache(["WURCS", "glycoct"]),
            )
        self.assertFalse(pd.isna(sugar_smiles[CHITOBIOSE_WURCS]))
        mock_glycoct.assert_not_called()

    def test_offline_failure_falls_back_to_live_chain(self):
        bad_wurcs = "not a wurcs string"
        with mock.patch.object(pac, "get_glycoct_from_wurcs", return_value="fake_glycoct") as mock_glycoct, \
             mock.patch.object(pac, "get_csdb_from_glycoct", return_value="fake_csdb") as mock_csdb, \
             mock.patch.object(pac, "get_smiles_from_csdb", return_value="C") as mock_smiles:
            sugar_smiles, *_ = pac.get_sugar_smiles_from_wurcs(
                [bad_wurcs],
                self._empty_cache(["glycoct", "csdb"]),
                self._empty_cache(["csdb", "descriptor"]),
                self._empty_cache(["WURCS", "glycoct"]),
            )
        mock_glycoct.assert_called_once()
        mock_csdb.assert_called_once()
        mock_smiles.assert_called_once()
        self.assertEqual(sugar_smiles[bad_wurcs], "C")


if __name__ == "__main__":
    unittest.main()
