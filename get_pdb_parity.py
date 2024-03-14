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

from utils import pdbe_sanitise_smiles
from rdkit.Chem.Draw import MolsMatrixToGridImage

def parity_score_smiles(row, threshold, progress, task, image_path): #, 
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
        scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "compound": None, "score" : 0, "error" : f"No biological compounds found for ligand", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0, parity_match: None}
        return scores_dict
    else:
        with redirect_stderr(f):
            try:
                #repeat canonicalisation to ensure best possible parity score
                ligand_rdkit, ligand_sanitisation = pdbe_sanitise_smiles(smiles, return_mol = True, return_sanitisation=True)
            except Exception as e:
                scores_dict = {"ec" : ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "cognate_ligand": cognate_ligand_id, "score" : 0, "error" : f"PDB Ligand error: {str(e)}", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0, parity_match: None}
                return scores_dict

        out = f.getvalue()
        
        ec = ",".join(row["entry"])
        try:
            rdkit_compound, compound_sanitisation = pdbe_sanitise_smiles(cognate_ligand_smiles, return_mol = True, return_sanitisation=True)
            if rdkit_compound in [None, np.nan]:
                scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "cognate_ligand": cognate_ligand_id, "score" : 0, "error" : f"RDKit compound not found for compound", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0, parity_match: None}
            else:
                if ligand_sanitisation and compound_sanitisation:
                    #repeat canonicalisation to ensure best possible parity score
                    parity = compare_molecules(ligand_rdkit, rdkit_compound, thresh = threshold)
                    #re-run the mcs to save the graph - or can we highlight the resulting atoms from the dictioanry returned by compare_molecules
                    score = parity.similarity_score
                    parity_match = parity.mapping
                    try:
                        mcs_graph = rdFMCS.FindMCS(
                            [mol1, mol2],
                            bondCompare=rdFMCS.BondCompare.CompareAny,
                            atomCompare=rdFMCS.AtomCompare.CompareAny,
                            timeout=40,
                            completeRingsOnly=True,
                        )
                        mcsMolsMatrix = [[mol1, mcs, mol2]]
                        mcsLegends = [[f'PDB Ligand UID: {int(pdb_ligand_id)}', f"Parity Score: {round(parity_score,3)}", f"Cognate Ligand UID: {int(cognate_ligand_id)}"]]
                        mcsHighlightAtomLists = [[(parity_match.keys()), (), (parity_match.values())]]
                        dos = Draw.MolDrawOptions()
                        dos.legendFraction = 0.10
                        img = Draw.MolsMatrixToGridImage(molsMatrix=mcsMolsMatrix, highlightAtomListsMatrix=mcsHighlightAtomLists, legendsMatrix=mcsLegends, returnPNG=True, drawOptions = dos)
                        with open(f'{image_path}/{int(pdb_ligand_id)}_{int(cognate_ligand_id)}_parity.png', 'wb') as f:
                            f.write(img.data)
                    except Exception: 
                        pass
                    
                    mol1_atom_count = ligand_rdkit.GetNumAtoms()
                    mol2_atom_count = rdkit_compound.GetNumAtoms()
                    matching_atoms = len(parity.mapping)
                    pdbl_subparity = matching_atoms/mol1_atom_count
                    bl_subparity = matching_atoms/mol2_atom_count
                    
                    cancelled = None
                    scores_dict = {"ec": ec, "pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "cognate_ligand": cognate_ligand_id, "score" : score, "error" : None, "cancelled" : cancelled, "pdbl_subparity": pdbl_subparity, "bl_subparity": bl_subparity, parity_match: parity_match}
                else:
                    scores_dict = {"ec": ec,"pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "cognate_ligand": cognate_ligand_id, "score" : 0, "error" : f"Sanitisation error", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0, parity_match: None}
        except Exception as e:
            scores_dict = {"ec": ec,"pdb_ligand" : pdb_ligand_id, "pdb_ligand_name": bl_name, "pdb_ligand_description": ligand_description, "cognate_ligand": cognate_ligand_id, "score" : 0, "error" : f"Parity error: {str(e)}", "cancelled" : None, "pdbl_subparity": 0, "bl_subparity": 0, parity_match: None}
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
image_path = f"{args.outdir}/parity_images"
Path(pickle_filename).parent.mkdir(parents=True, exist_ok=True)
Path(image_path).mkdir(parents=True, exist_ok=True)
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
                futures.append(executor.submit(parity_score_smiles, row, args.threshold, progress,task, image_path))
            for future in concurrent.futures.as_completed(futures):
                results.extend([future.result()])
                futures.remove(future)

        results_df = pd.DataFrame(results)
        # Save the smiles_ec_pairs for the chunk as a pickle file
        results_df.to_pickle(pickle_filename)
end = time.time()
print(f"Time taken: {end - start} seconds")