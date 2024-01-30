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
#need to track whether a molecules MCS was timeouted or not, and then we can collect and repeat these structures 

def get_parity_score(mol1, mol2, print_structures = False, ringmatches = False,returnmcs = False, timeout = 300):
    # Get the MCS of the two molecules
    mcs_result = rdFMCS.FindMCS(
        [mol1, mol2], 
        matchValences=False,
        ringMatchesRingOnly= ringmatches,
        atomCompare=rdFMCS.AtomCompare.CompareAny, 
        bondCompare=rdFMCS.BondCompare.CompareAny, timeout = timeout)
    if mcs_result.canceled:
        cancelled = True
    else:
        cancelled = False
    # Convert the SMARTS string from the MCS result into a molecule
    mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)

    # Map the MCS back onto the original molecules and get the atom indices
    mol1_match = mol1.GetSubstructMatch(mcs_mol)
    mol2_match = mol2.GetSubstructMatch(mcs_mol)
    
    if print_structures:
        # Draw the molecules and the MCS
        Draw.MolToMPL(mol1, size=(200,200), kekulize=True, title='Molecule 1')
        Draw.MolToMPL(mol2, size=(200,200), kekulize=True, title='Molecule 2')
        Draw.MolToMPL(mcs_mol, size=(200,200), kekulize=True, title='MCS')
    # Compare the atom types at each matched position
    Nsim = 0
    for idx1, idx2 in zip(mol1_match, mol2_match):
        atom1 = mol1.GetAtomWithIdx(idx1)
        atom2 = mol2.GetAtomWithIdx(idx2)
        if atom1.GetAtomicNum() == atom2.GetAtomicNum():
            Nsim += 1

    # Compute PARITY similarity score
    NB = mol1.GetNumAtoms()
    NC = mol2.GetNumAtoms()
    score = Nsim / (NB + NC - Nsim)
    if returnmcs:
        return score, cancelled, mcs_mol
    else:
        return score, cancelled

    
import io
from contextlib import redirect_stderr

def parity_score_inchi(pdb_ligand_id, inchi, ec, bl_name, compound_df, timeout = 300):
    f = io.StringIO()
    scores = []
    
    compound_df_subset = compound_df.loc[(compound_df.entry.isin(ec)), ["entry", "canonical_smiles", "ROMol"]]
    compound_df_subset = compound_df_subset.groupby("canonical_smiles").agg({"entry" : list, "ROMol" : "first"}).reset_index()
    ec = ",".join(ec)
    if len(compound_df_subset) == 0:
        score_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "compound": None, "score" : 0, "error" : f"No biological compounds found for ligand", "cancelled" : None}
        scores.append(score_dict)
    else:
        with redirect_stderr(f):
            try:
                ligand_rdkit = Chem.inchi.MolFromInchi(inchi)
            except Exception as e:
                for compound in compound_list:
                    score_dict = {"ec" : ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "compound": compound, "score" : 0, "error" : f"PDB Ligand error: {str(e)}", "cancelled" : None}
                    scores.append(score_dict)
                return scores

        out = f.getvalue()
        for idx, row in compound_df_subset.iterrows():
            compound = row["canonical_smiles"]
            ec = ",".join(row["entry"])
            try:
                rdkit_compound = row["ROMol"]
                if rdkit_compound in [None, np.nan]:
                    scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "compound": compound, "score" : 0, "error" : f"RDKit compound not found for compound", "cancelled" : None}
                else:
                    score, cancelled = get_parity_score(ligand_rdkit, rdkit_compound, returnmcs=False, timeout = timeout)
                    scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "compound": compound, "score" : score, "error" : None, "cancelled" : cancelled}
                scores.append(scores_dict)
            except Exception as e:
                scores_dict = {"ec": ec,"pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "compound": compound, "score" : 0, "error" : f"Parity error: {str(e)}", "cancelled" : None}
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
parser.add_argument('--chunk', metavar = '', type = int, 
    help = "")
parser.add_argument('--chunk_size', metavar = '', type = int, default= 100,
    help = "")

args = parser.parse_args()

all_chem_descriptors_ligands_unique_pairs = pd.read_pickle(f"{args.processed_ligands_file}")

cognate_ligands_df = pd.read_pickle(f"{args.cognate_ligands_file}")
chunk_index = args.chunk
chunk_size = args.chunk_size


start_index = chunk_index * chunk_size
end_index = start_index + chunk_size
chunk_df = all_chem_descriptors_ligands_unique_pairs[start_index:end_index]
pickle_filename = f"{args.outdir}/parity_calcs_chunk_{chunk_index}.pkl"


if not os.path.exists(pickle_filename):
    chunk_results = []
    for index, row in chunk_df.iterrows():
        bl_name = row['bl_name']
        x = row['ligand_entity_id']
        y = row['descriptor']
        z = row['ec_list']
        scores_list = parity_score_inchi(x, y, z, bl_name, cognate_ligands_df, timeout = 30)
        chunk_results.extend(scores_list)
        print(f"Chunk {args.chunk}, row {index}")
    # Save the inchi_ec_pairs for the chunk as a pickle file
    with open(pickle_filename, "wb") as file:
        pickle.dump(chunk_results, file)
