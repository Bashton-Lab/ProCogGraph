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

import argparse
import snakemake

parser = argparse.ArgumentParser(description = 'TO DO')
parser.add_argument('--pdb_ligands_file', default = 'cath_scop_ligands_to_score.pkl', type = str, 
    help = "")
parser.add_argument('--cognate_ligands_file', default = 'final_kegg_compound_reaction_enzyme_df.pkl', type = str, 
    help = "")
parser.add_argument('--all_pdb_ligands', default = 'chains_domains_chem_descriptors_cath_scop.pkl', type = str, 
    help = "")
parser.add_argument('--outdir', type = str, default = "parity_calcs", 
    help = ""),
parser.add_argument('--chunk_size', type = int, default= 1000, 
    help = "")
parser.add_argument('--threads', type = int, default= 20, 
    help = "")
parser.add_argument('--snakefile', default = 'parity.smk', type = str,
    help = "")

args = parser.parse_args()

chains_domains_chem_descriptors = pd.read_pickle(args.pdb_ligands_file)
        
num_chunks = math.ceil(len(chains_domains_chem_descriptors) / args.chunk_size)
config = { 
    "outdir" : args.outdir,
    "num_chunks" : num_chunks,
    "threads" : args.threads,
    "chunk_size" : args.chunk_size,
    "prepared_ligands" : args.pdb_ligands_file,
    "all_pdb_ligands" : args.all_pdb_ligands,
    "cognate_ligands_file": args.cognate_ligands_file}

status = snakemake.snakemake(
    args.snakefile, 
    printshellcmds=False,
    config=config,
    cores = args.threads)