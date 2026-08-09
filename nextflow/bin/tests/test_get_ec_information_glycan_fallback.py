#!/usr/bin/env python

"""
get_ec_information.py's Context A glycan-resolution block is inline
script code (not a function), so it can't be called directly in a test.
This instead replicates its exact offline-fallback logic
(missing_smiles_mask / .apply(get_smiles_from_wurcs_offline)) against a
small synthetic dataframe, to validate that pattern in isolation without
needing the full pipeline's upstream KEGG/Rhea/GlyTouCan machinery.

    python3 nextflow/bin/tests/test_get_ec_information_glycan_fallback.py
"""

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils import get_smiles_from_wurcs_offline

CHITOBIOSE_WURCS = "WURCS=2.0/1,2,1/[a2122h-1b_1-5_2*NCC/3=O]/1-1/a4-b1"


class TestOfflineFallbackMasking(unittest.TestCase):

    def test_only_missing_rows_with_a_wurcs_value_are_backfilled(self):
        df = pd.DataFrame({
            "compound_id": ["G1", "G2", "G3", "G4"],
            # G1: live chain already resolved it -> must not be overwritten
            # G2: live chain failed, has a wurcs -> should be backfilled offline
            # G3: live chain failed, no wurcs at all -> stays nan
            # G4: live chain failed, wurcs is unparseable -> stays nan
            "smiles": ["C(already resolved)", np.nan, np.nan, np.nan],
            "wurcs": [CHITOBIOSE_WURCS, CHITOBIOSE_WURCS, np.nan, "not a wurcs string"],
        })

        # exact logic from get_ec_information.py
        missing_smiles_mask = df["smiles"].isna() & df["wurcs"].notna()
        df.loc[missing_smiles_mask, "smiles"] = df.loc[missing_smiles_mask, "wurcs"].apply(get_smiles_from_wurcs_offline)

        self.assertEqual(df.loc[0, "smiles"], "C(already resolved)")
        self.assertIsInstance(df.loc[1, "smiles"], str)
        self.assertNotEqual(df.loc[1, "smiles"], "")
        self.assertTrue(pd.isna(df.loc[2, "smiles"]))
        self.assertTrue(pd.isna(df.loc[3, "smiles"]))


if __name__ == "__main__":
    unittest.main()
