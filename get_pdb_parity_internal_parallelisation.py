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

def parity_score_smiles(pdb_ligand_id, smiles, ec, bl_name, ligand_description, compound_df_subset, threshold):
    f = io.StringIO()
    scores = []
    
    ec = ",".join(ec)
    if len(compound_df_subset) == 0:
        score_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": None, "score" : 0, "error" : f"No biological compounds found for ligand", "cancelled" : None}
        scores.append(score_dict)
    else:
        with redirect_stderr(f):
            try:
                #repeat canonicalisation to ensure best possible parity score
                ligand_rdkit = Chem.MolFromSmiles(Chem.CanonSmiles(smiles))
            except Exception as e:
                for idx, row in compound_df_subset.iterrows():
                    compound = row["canonical_smiles"]
                    score_dict = {"ec" : ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": compound, "score" : 0, "error" : f"PDB Ligand error: {str(e)}", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0}
                    scores.append(score_dict)
                return scores

        out = f.getvalue()
        for idx, row in compound_df_subset.iterrows():
            compound = row["canonical_smiles"]
            ec = ",".join(row["entry"])
            try:
                rdkit_compound = row["ROMol"]
                if rdkit_compound in [None, np.nan]:
                    scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": compound, "score" : 0, "error" : f"RDKit compound not found for compound", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0}
                else:
                    #repeat canonicalisation to ensure best possible parity score
                    rdkit_compound = Chem.MolFromSmiles(Chem.CanonSmiles(row["canonical_smiles"]))
                    parity = compare_molecules(ligand_rdkit, rdkit_compound, thresh = threshold)
                    score = parity.similarity_score
                    mol1_atom_count = ligand_rdkit.GetNumAtoms()
                    mol2_atom_count = rdkit_compound.GetNumAtoms()
                    matching_atoms = len(parity.mapping)
                    pdbl_subparity = matching_atoms/mol1_atom_count
                    bl_subparity = matching_atoms/mol2_atom_count
                    cancelled = None
                    #score, pdbl_subparity, bl_subparity, cancelled = get_parity_score(ligand_rdkit, rdkit_compound, returnmcs=False, timeout = timeout)
                    scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": compound, "score" : score, "error" : None, "cancelled" : cancelled, "pdbl_subparity": pdbl_subparity, "bl_subparity": bl_subparity}
                scores.append(scores_dict)
            except Exception as e:
                scores_dict = {"ec": ec,"pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": compound, "score" : 0, "error" : f"Parity error: {str(e)}", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0}
                scores.append(scores_dict)
    
    return scores


import argparse

parser = argparse.ArgumentParser(description = 'TO DO')
parser.add_argument('--processed_ligands_file', metavar = '', type = str, 
    help = "")
parser.add_argument('--cognate_ligands_file', metavar = '', type = str, 
    help = "")
parser.add_argument('--outdir', metavar = 'parity_calcs', type = str, 
    help = "")

args = parser.parse_args()

all_chem_descriptors_ligands_unique_pairs = pd.read_pickle(f"{args.processed_ligands_file}")
all_chem_descriptors_ligands_unique_pairs = all_chem_descriptors_ligands_unique_pairs.head(1000)
cognate_ligands_df = pd.read_pickle(f"{args.cognate_ligands_file}")
chunk_index = args.chunk
chunk_size = args.chunk_size

pickle_filename = f"{args.outdir}/all_parity_calcs.pkl"

Path(pickle_filename).parent.mkdir(parents=True, exist_ok=True)
if not os.path.exists(pickle_filename):
    results = []
    with Progress() as progress:
        task = progress.add_task("[red]Calculating parity scores...", total=len(all_chem_descriptors_ligands_unique_pairs))
        
        with concurrent.futures.ThreadPoolExecutor(max_workers = 50) as executor:
            futures = []
            for index, row in all_chem_descriptors_ligands_unique_pairs.iterrows():
                cognate_ligands_df_subset = cognate_ligands_df.loc[(cognate_ligands_df.entry.isin(row.ec_list)), ["entry", "canonical_smiles", "ROMol"]].copy()
                cognate_ligands_df_subset = cognate_ligands_df_subset.groupby("canonical_smiles").agg({"entry" : list, "ROMol" : "first"}).reset_index()
                futures.append(executor.submit(parity_score_smiles, row, cognate_ligands_df_subset, progress,task) for _, row in all_chem_descriptors_ligands_unique_pairs.iterrows())
            for future in concurrent.futures.as_completed(futures):
                results.extend(future.result())
                del futures[future]

        results_df = pd.DataFrame(results)
        # Save the smiles_ec_pairs for the chunk as a pickle file
        results_df.to_pickle(pickle_filename)