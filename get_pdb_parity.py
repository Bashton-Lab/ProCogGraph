#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import requests
from Bio.KEGG import Enzyme
import matplotlib
import io
import time
import pickle
import re
from rich.progress import Progress
from rdkit import Chem
from collections import Counter
import rdkit.Chem as Chem
from rdkit.Chem import BRICS,Recap, rdFMCS
from rdkit.Chem import Draw

import pickle
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

def parity_score_smiles(row, threshold, progress, task): #, 
    ec = row.entry
    pdb_ligand_id = row.pdb_ligand_id
    smiles = row.smiles
    bl_name = row.bl_name
    ligand_description = row.ligand_description
    cognate_ligand_smiles = row.canonical_smiles
    f = io.StringIO()
    
    ec = ",".join(ec)
    if cognate_ligand_smiles == None:
        scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": None, "score" : 0, "error" : f"No biological compounds found for ligand", "cancelled" : None}
        return scores_dict
    else:
        with redirect_stderr(f):
            try:
                #repeat canonicalisation to ensure best possible parity score
                ligand_rdkit = Chem.MolFromSmiles(Chem.CanonSmiles(smiles))
            except Exception as e:
                scores_dict = {"ec" : ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": cognate_ligand_smiles, "score" : 0, "error" : f"PDB Ligand error: {str(e)}", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0}
                return scores_dict

        out = f.getvalue()
        
        ec = ",".join(row["entry"])
        try:
            rdkit_compound = Chem.MolFromSmiles(Chem.CanonSmiles(cognate_ligand_smiles))
            if rdkit_compound in [None, np.nan]:
                scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": cognate_ligand_smiles, "score" : 0, "error" : f"RDKit compound not found for compound", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0}
            else:
                #repeat canonicalisation to ensure best possible parity score
                parity = compare_molecules(ligand_rdkit, rdkit_compound, thresh = threshold)
                score = parity.similarity_score
                mol1_atom_count = ligand_rdkit.GetNumAtoms()
                mol2_atom_count = rdkit_compound.GetNumAtoms()
                matching_atoms = len(parity.mapping)
                pdbl_subparity = matching_atoms/mol1_atom_count
                bl_subparity = matching_atoms/mol2_atom_count
                cancelled = None
                #score, pdbl_subparity, bl_subparity, cancelled = get_parity_score(ligand_rdkit, rdkit_compound, returnmcs=False, timeout = timeout)
                scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": cognate_ligand_smiles, "score" : score, "error" : None, "cancelled" : cancelled, "pdbl_subparity": pdbl_subparity, "bl_subparity": bl_subparity}
        except Exception as e:
            scores_dict = {"ec": ec,"pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": cognate_ligand_smiles, "score" : 0, "error" : f"Parity error: {str(e)}", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0}
    progress.update(task, advance=1)
    return scores_dict

def get_compound_pairs(row, cognate_ligands_df):
    ec = row.ec_list
    ec = [ec_number for sublist in (ec_list.split(',') for ec_list in ec) for ec_number in sublist] #split each list of ECs into individual ECs
    pdb_ligand_id = row.ligand_entity_id
    smiles = row.descriptor
    bl_name = row.bl_name
    ligand_description = row.description
    cognate_ligands_df_subset = cognate_ligands_df.loc[(cognate_ligands_df.entry.isin(ec)), ["entry", "canonical_smiles", "ROMol"]].copy()
    cognate_ligands_df_subset = cognate_ligands_df_subset.groupby("canonical_smiles").agg({"entry" : list, "ROMol" : "first"}).reset_index()
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
parser.add_argument('--outdir', metavar = 'parity_calcs', type = str, 
    help = "")
parser.add_argument('--threshold', metavar = '', type = float, default = 0.1,
    help = "")
parser.add_argument('--threads', metavar = '', type = int, default = 1,
    help = "")

start = time.time()
args = parser.parse_args()
RDLogger.DisableLog('rdApp.*') 
all_chem_descriptors_ligands_unique_pairs = pd.read_pickle(f"{args.processed_ligands_file}")
cognate_ligands_df = pd.read_pickle(f"{args.cognate_ligands_file}")

pickle_filename = f"{args.outdir}/all_parity_calcs.pkl"

Path(pickle_filename).parent.mkdir(parents=True, exist_ok=True)
if not os.path.exists(pickle_filename):
    all_pairs = []
    for index, row in all_chem_descriptors_ligands_unique_pairs.iterrows():
        pairs = get_compound_pairs(row, cognate_ligands_df)
        all_pairs.append(pairs)

    all_pairs_df = pd.concat(all_pairs, ignore_index = True)
    results = []
    with Progress() as progress:
        task = progress.add_task("[red]Calculating parity scores...", total=len(all_pairs_df))
        
        with concurrent.futures.ThreadPoolExecutor(args.threads) as executor:
            futures = []
            for index, row in all_pairs_df.iterrows():
                futures.append(executor.submit(parity_score_smiles, row, args.threshold, progress,task))
            for future in concurrent.futures.as_completed(futures):
                results.extend([future.result()])
                futures.remove(future)

        results_df = pd.DataFrame(results)
        # Save the smiles_ec_pairs for the chunk as a pickle file
        results_df.to_pickle(pickle_filename)
end = time.time()
print(f"Time taken: {end - start} seconds")