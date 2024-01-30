import math
import os
import pandas as pd
import pickle
from rich.progress import Progress
from functools import partial
#python3 snakemake_ligands_df.py --pdb_ligands_file /raid/MattC/repos/CognateLigandProject/pdbe_graph_files/bound_entities_to_score.pkl --cognate_ligands /raid/MattC/repos/CognateLigandProject/biological_ligands/biological_ligands_df.pkl --outdir bound_entities_parity_3 --chunk_size 10 --threads 70 --snakefile parity.smk

rule all:
    input:
        config["outdir"] + "/all_parity_calcs.pkl"

rule calculate_parity:
    input:
        infile = config["prepared_ligands"],
        cognate_ligands = config["cognate_ligands_file"]
    output:
        chunk_results= config["outdir"] + "/parity_calcs_chunk_{chunk_index}.pkl"
    params:
        chunk_size=config["chunk_size"],
        outdir = config["outdir"]
    shell:
        "python3 get_pdb_parity.py --processed_ligands_file {input.infile} --cognate_ligands_file {input.cognate_ligands} --chunk {wildcards.chunk_index} --chunk_size {params.chunk_size} --outdir {params.outdir}"

rule merge_results:
    input:
        chunk_results=expand(config['outdir'] + "/parity_calcs_chunk_{chunk_index}.pkl", chunk_index=range(config["num_chunks"]))
    output:
        config["outdir"] + "/all_parity_calcs.pkl"
    shell:
        "python3 merge_parity_chunks.py --chunk_list {input.chunk_results} --output_pkl {output} --ligands_score_subset {config[prepared_ligands]}"