#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import requests
from Bio.KEGG import Enzyme
import io
import time
import pickle
import re
from rdkit import Chem
from collections import Counter

from rdkit.Chem import rdFMCS
import os
import math
import numpy as np
from math import isnan
import concurrent.futures
from pathlib import Path
#need to track whether a molecules MCS was timeouted or not, and then we can collect and repeat these structures 

from pdbeccdutils.computations.parity_method import compare_molecules
    
import io
from contextlib import redirect_stderr

from utils import pdbe_sanitise_smiles
from rdkit.Chem.Draw import MolsMatrixToGridImage

def parity_score_smiles(row, threshold): #, 
    ec = row.entry
    pdb_ligand_id = row.pdb_ligand_id
    smiles = row.smiles
    bl_name = row.bl_name
    ligand_description = row.ligand_description
    cognate_ligand_smiles = row.canonical_smiles
    cognate_ligand_id = row.uniqueID
    f = io.StringIO()
    
    ec = ",".join(ec)

    if cognate_ligand_smiles == None:
        scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "pdb_ligand_smiles": smiles, "cognate_ligand": cognate_ligand_id, "cognate_ligand_smiles": None, "score" : 0, "error" : f"No biological compounds found for ligand", "pdbl_subparity": 0, "bl_subparity": 0, "parity_match": None, "parity_smarts": None, "threshold": threshold}
        return scores_dict
    else:
        with redirect_stderr(f):
            try:
                #repeat canonicalisation to ensure best possible parity score
                ligand_rdkit, ligand_sanitisation = pdbe_sanitise_smiles(smiles, return_mol = True, return_sanitisation=True)
            except Exception as e:
                scores_dict = {"ec" : ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "pdb_ligand_smiles": smiles, "cognate_ligand": cognate_ligand_id, "cognate_ligand_smiles": None, "score" : 0, "error" : f"PDB Ligand error: {str(e)}", "pdbl_subparity": 0, "bl_subparity": 0, "parity_match": None, "parity_smarts": None, "threshold": threshold}
                return scores_dict

        out = f.getvalue()
        
        ec = ",".join(row["entry"])
        try:
            rdkit_compound, compound_sanitisation = pdbe_sanitise_smiles(cognate_ligand_smiles, return_mol = True, return_sanitisation=True)
            if rdkit_compound in [None, np.nan]:
                scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "pdb_ligand_smiles": smiles, "cognate_ligand": cognate_ligand_id, "cognate_ligand_smiles": None, "score" : 0, "error" : f"RDKit compound not found for compound", "pdbl_subparity": 0, "bl_subparity": 0, "parity_match": None, "parity_smarts": None, "threshold": threshold}
            else:
                if ligand_sanitisation and compound_sanitisation:
                    #repeat canonicalisation to ensure best possible parity score
                    parity = compare_molecules(ligand_rdkit, rdkit_compound, thresh = threshold)
                    #re-run the mcs to save the graph - or can we highlight the resulting atoms from the dictioanry returned by compare_molecules
                    score = parity.similarity_score
                    parity_match = parity.mapping
                    mcs_graph = rdFMCS.FindMCS(
                                [ligand_rdkit, rdkit_compound],
                                bondCompare=rdFMCS.BondCompare.CompareAny,
                                atomCompare=rdFMCS.AtomCompare.CompareAny,
                                timeout=40,
                                completeRingsOnly=True,
                            )
                    mcs_smarts = mcs_graph.smartsString
                    
                    mol1_atom_count = ligand_rdkit.GetNumAtoms()
                    mol2_atom_count = rdkit_compound.GetNumAtoms()
                    matching_atoms = len(parity.mapping)
                    pdbl_subparity = matching_atoms/mol1_atom_count
                    bl_subparity = matching_atoms/mol2_atom_count
                    
                    scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "pdb_ligand_smiles": smiles, "cognate_ligand": cognate_ligand_id, "cognate_ligand_smiles": cognate_ligand_smiles, "score" : score, "error" : None, "pdbl_subparity": pdbl_subparity, "bl_subparity": bl_subparity, "parity_match": parity_match, "parity_smarts": mcs_smarts, "threshold": threshold}
                else:
                    scores_dict = {"ec": ec,"pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "pdb_ligand_smiles": smiles,  "cognate_ligand": cognate_ligand_id, "cognate_ligand_smiles": cognate_ligand_smiles, "score" : 0, "error" : f"Sanitisation error", "pdbl_subparity": 0, "bl_subparity": 0, "parity_match": None, "parity_smarts": None, "threshold": threshold}
        except Exception as e:
            scores_dict = {"ec": ec,"pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "pdb_ligand_smiles": smiles, "cognate_ligand": cognate_ligand_id, "cognate_ligand_smiles": cognate_ligand_smiles, "score" : 0, "error" : f"Parity error: {str(e)}", "pdbl_subparity": 0, "bl_subparity": 0, "parity_match": None, "parity_smarts": None, "threshold": threshold}
    return scores_dict

def get_compound_pairs(row, cognate_ligands_df):
    ec = row.ec_list
    ec = [ec_number for sublist in (ec_list.split(',') for ec_list in ec) for ec_number in sublist] #split each list of ECs into individual ECs
    pdb_ligand_id = row.ligand_entity_id
    smiles = row.descriptor
    bl_name = row.bl_name
    ligand_description = row.description
    cognate_ligands_df_subset = cognate_ligands_df.loc[(cognate_ligands_df.entry.isin(ec)), ["entry", "canonical_smiles", "uniqueID"]].copy()
    cognate_ligands_df_subset = cognate_ligands_df_subset.groupby("canonical_smiles").agg({"entry" : list, "uniqueID": "first"}).reset_index()
    cognate_ligands_df_subset[["pdb_ligand_id", "smiles", "bl_name", "ligand_description"]] = [pdb_ligand_id, smiles, bl_name, ligand_description]
    return cognate_ligands_df_subset

import argparse
from rdkit import RDLogger
import time

parser = argparse.ArgumentParser(description = 'TO DO')
parser.add_argument('--processed_ligands_file', metavar = '', type = str, 
    help = "")
parser.add_argument('--cognate_ligands_file', metavar = '', type = str, 
    help = "")
parser.add_argument('--threshold', metavar = '', type = float, default = 0.1,
    help = "")
parser.add_argument('--threads', metavar = '', type = int, default = 1,
    help = "")
parser.add_argument('--cache', metavar = '', type = str,
    help = "")

start = time.time()
args = parser.parse_args()
RDLogger.DisableLog('rdApp.*') 
all_chem_descriptors_ligands_unique_pairs = pd.read_pickle(f"{args.processed_ligands_file}")
cognate_ligands_df = pd.read_pickle(f"{args.cognate_ligands_file}")

pickle_filename = f"all_parity_calcs.pkl"

if not os.path.exists(pickle_filename):
    all_pairs = []
    for index, row in all_chem_descriptors_ligands_unique_pairs.iterrows():
        pairs = get_compound_pairs(row, cognate_ligands_df)
        all_pairs.append(pairs)

    all_pairs_df = pd.concat(all_pairs, ignore_index = True)
    if args.cache:
        cache_df = pd.read_pickle(f"{args.cache}")
        cache_df = cache_df.drop_duplicates(subset = ["pdb_ligand_smiles", "cognate_ligand_smiles"]) #double in case not done in new cache results below (left over from first run)
    else:
        cache_df = None

    all_pairs_df2 = all_pairs_df.merge(cache_df, left_on = ["smiles", "canonical_smiles"], right_on = ["pdb_ligand_smiles", "cognate_ligand_smiles"], how = "left", validate = "many_to_one", indicator = True)
    pre_calculated = all_pairs_df2.loc[all_pairs_df2._merge == "both", ['entry', 'pdb_ligand_id', 'bl_name', 'ligand_description', 'uniqueID', 'score', 'error', 'pdbl_subparity', 'bl_subparity', 'parity_match', 'parity_smarts', 'pdb_ligand_smiles', 'cognate_ligand_smiles', 'threshold']].reset_index(drop = True)
    pre_calculated.rename(columns = {"entry":"ec", "ligand_description": "pdb_ligand_description", "uniqueID": "cognate_ligand", "pdb_ligand_id": "pdb_ligand", "bl_name": "pdb_ligand_name"}, inplace = True)
    pre_calculated["ec"] = pre_calculated.ec.str.join(",")
    pre_calculated["cache"] = True
    to_calculate = all_pairs_df2.loc[all_pairs_df2._merge == "left_only"]
    print(f"Pre-calculated: {len(pre_calculated)}")
    print(f"To calculate: {len(to_calculate)}")
    results = []
    if len(to_calculate) > 0:
        with concurrent.futures.ThreadPoolExecutor(args.threads) as executor:
            futures = []
            for index, row in to_calculate.iterrows():
                futures.append(executor.submit(parity_score_smiles, row, args.threshold))
            for future in concurrent.futures.as_completed(futures):
                results.extend([future.result()])
                futures.remove(future)

        results_df = pd.DataFrame(results)
        results_df = pd.concat([results_df, pre_calculated], ignore_index = True).reset_index(drop = True)
        # Save the smiles_ec_pairs for the chunk as a pickle file
        results_df.to_pickle(pickle_filename)

        #save the cache scores
        cache_results = results_df[["pdb_ligand_smiles", "cognate_ligand_smiles", "threshold", "score", "error", "pdbl_subparity", "bl_subparity", "parity_match", "parity_smarts"]].drop_duplicates(subset = ["pdb_ligand_smiles", "cognate_ligand_smiles"])
        #when dropping duplicates here consider that it needs to drop the max theshold duplicate - better to do this in a groupby?
        cache_results = pd.concat([cache_df, cache_results], ignore_index = True).drop_duplicates(subset = ["pdb_ligand_smiles", "cognate_ligand_smiles", "threshold"])
        cache_results.to_pickle(f"cache_parity_calcs_new.pkl")
    else:
        results_df = pre_calculated
        # Save the smiles_ec_pairs for the chunk as a pickle file
        results_df.to_pickle(pickle_filename)
        cache_df.to_pickle(f"cache_parity_calcs_new.pkl") #output the original cache as new if no new calculations are done


end = time.time()
print(f"Time taken: {end - start} seconds")