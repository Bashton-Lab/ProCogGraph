import matplotlib.pyplot as plt
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
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS
import pubchempy as pcp
from bs4 import BeautifulSoup
from urllib.parse import quote
import json
from pathlib import Path
import math
import os

from ast import literal_eval
import argparse

def assign_ownership_percentile_categories(ligands_df, unique_id = "uniqueID", domain_grouping_key = "cath_domain"):
    ligands_df["total_contact_counts"] = ligands_df.groupby([unique_id]).transform("size")
    ligands_df["domain_contact_counts"] = ligands_df.groupby([unique_id, domain_grouping_key]).transform("size")
    hbond_counts = ligands_df.explode('contact_type').groupby([unique_id, domain_grouping_key])['contact_type'].apply(lambda x: (x == 'hbond').sum()).rename("domain_hbond_counts").reset_index()
    ligands_df = ligands_df.merge(hbond_counts, how = "left", on = [unique_id, domain_grouping_key], indicator = True)
    assert(len(ligands_df.loc[ligands_df._merge != "both"]) == 0)
    ligands_df.drop(columns = ["_merge"], inplace = True)
    covalent_counts = ligands_df.explode('contact_type').groupby([unique_id, domain_grouping_key])['contact_type'].apply(lambda x: (x == 'covalent').sum()).rename("domain_covalent_counts").reset_index()
    ligands_df = ligands_df.merge(covalent_counts, how = "left", on = [unique_id, domain_grouping_key], indicator = True)
    assert(len(ligands_df.loc[ligands_df._merge != "both"]) == 0)
    ligands_df.drop(columns = ["_merge"], inplace = True)
    ligands_df["domain_hbond_perc"] = ligands_df.domain_hbond_counts / ligands_df.total_contact_counts
    ligands_df["domain_contact_perc"] = ligands_df.domain_contact_counts / ligands_df.total_contact_counts
    ligands_df["domain_ownership"] = np.where(
        ligands_df["domain_covalent_counts"] > 0, "covalent",
        np.where(
            ligands_df["domain_contact_perc"] == 1, "unique",
            np.where(
                ligands_df["domain_contact_perc"] >= 0.7, "dominant",
                np.where(
                    (ligands_df["domain_contact_perc"] >= 0.3)
                    & (ligands_df["domain_contact_perc"] < 0.7), "partner",
                    np.where(
                        ligands_df["domain_contact_perc"] < 0.3, "minor", np.nan)
                )
            )
        )
    )
    
    return ligands_df

def sorted_set(x):
    return sorted(set(x))

def main():
    parser = argparse.ArgumentParser(description = 'TO DO')
    parser.add_argument('--cath_bl_residue_interactions_file', metavar = '', type = str,
                        help = "cath_pdb_residue_interactions_distinct_bl_ec.csv.gz")
    parser.add_argument('--scop_bl_residue_interactions_file', metavar = '', type = str,
                        help = "scop_pdb_residue_interactions_distinct_bl_ec.csv.gz")
    parser.add_argument('--interpro_bl_residue_interactions_file', metavar = '', type = str,
                        help = "interpro_pdb_residue_interactions_distinct_bl_ec.csv.gz")
    parser.add_argument('--cath_sugar_residue_interactions_file', metavar = '', type = str,
                        help = "cath_pdb_residue_interactions_distinct_sugar_ec.csv.gz")
    parser.add_argument('--scop_sugar_residue_interactions_file', metavar = '', type = str,
                        help = "scop_pdb_residue_interactions_distinct_sugar_ec.csv.gz")
    parser.add_argument('--interpro_sugar_residue_interactions_file', metavar = '', type = str,
                        help = "interpro_pdb_residue_interactions_distinct_sugar_ec.csv.gz")
    parser.add_argument('--outdir', metavar = '', type = str,
                        help = "output directory")
    parser.add_argument('--scop_domains_info_file', metavar = '', type = str,
                        help = "dir.cla.scop.1_75.txt")
    parser.add_argument('--scop_descriptions_file', metavar = '', type = str,
                        help = "dir.des.scop.1_75.txt")

    args = parser.parse_args()
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    cath_bl_residue_df = pd.read_csv(f"{args.cath_bl_residue_interactions_file}", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    cath_bl_residue_df.rename(columns = {"bound_ligand_name": "name"}, inplace = True)
    cath_sugar_residue_df = pd.read_csv(f"{args.cath_sugar_residue_interactions_file}", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    cath_sugar_residue_df.rename(columns = {"ligand_entity_description": "name"}, inplace = True)
    
    cath_bl_residue_df["contact_type"] = cath_bl_residue_df["contact_type"].apply(literal_eval)
    cath_bl_residue_df_domains = assign_ownership_percentile_categories(cath_bl_residue_df.copy(), "uniqueID", "cath_domain")
    cath_cols = ["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords", "protein_entity_id", "protein_entity_ec", "ec_list", "cath_domain", "cath_class", "cath_architecture", "cath_topology","cath_homology", "cath_name", "cath_code", "bm_uniqids", "bound_molecule_display_id", "name", "uniqueID", "bound_ligand_id", "ligand_uniqueID", "type","total_contact_counts", "domain_contact_counts", "domain_hbond_counts", "domain_contact_perc", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "protein_chain_id", "bound_ligand_struct_asym_id"]
    cath_bl_residue_df_domains["pdb_descriptor"] = cath_bl_residue_df_domains["pdb_descriptor"].fillna("")
    cath_bl_domains = cath_bl_residue_df_domains.groupby(cath_cols, dropna = False).agg({"pdb_residue_auth_id": sorted_set, "bound_ligand_auth_id": sorted_set}).reset_index()
    cath_bl_domains["pdb_residue_auth_id"] = cath_bl_domains["pdb_residue_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    cath_bl_domains["bound_ligand_auth_id"] = cath_bl_domains["bound_ligand_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    cath_bl_domains.drop_duplicates(inplace = True)
    cath_bl_domains.to_csv(f"{args.outdir}/cath_bl_domain_ownership.csv", index = False)

    cath_sugar_residue_df["contact_type"] = cath_sugar_residue_df["contact_type"].apply(literal_eval)
    cath_sugar_residue_df_domains = assign_ownership_percentile_categories(cath_sugar_residue_df.copy(), "uniqueID", "cath_domain")
    cath_sugar_domains = cath_sugar_residue_df_domains.groupby(cath_cols, dropna = False).agg({"pdb_residue_auth_id": sorted_set, "bound_ligand_auth_id": sorted_set}).reset_index()
    cath_sugar_domains["pdb_residue_auth_id"] = cath_sugar_domains["pdb_residue_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    cath_sugar_domains["bound_ligand_auth_id"] = cath_sugar_domains["bound_ligand_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    cath_sugar_domains.drop_duplicates(inplace = True)
    cath_sugar_domains.to_csv(f"{args.outdir}/cath_sugar_domain_ownership.csv", index = False)

    cath_combined_domains = pd.concat([cath_bl_domains, cath_sugar_domains])
    cath_combined_domains.to_csv(f"{args.outdir}/cath_combined_domain_ownership.csv", index = False)


    scop_bl_residue_df = pd.read_csv(f"{args.scop_bl_residue_interactions_file}", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    scop_bl_residue_df.rename(columns = {"bound_ligand_name": "name"}, inplace = True)
    scop_bl_residue_df["contact_type"] = scop_bl_residue_df["contact_type"].apply(literal_eval)
    scop_bl_residue_df_domains = assign_ownership_percentile_categories(scop_bl_residue_df.copy(), "uniqueID", "scop_id")

    scop_cols = ["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords", "protein_entity_id", "protein_entity_ec", "ec_list", "scop_sunid","scop_description", "scop_sccs", "scop_class_id", "scop_fold_id", "scop_superfamily_id", "scop_id", "bm_uniqids", "bound_molecule_display_id", "name", "uniqueID", "bound_ligand_id", "ligand_uniqueID", "type","total_contact_counts", "domain_contact_counts", "domain_hbond_counts", "domain_contact_perc", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "protein_chain_id", "bound_ligand_struct_asym_id"]
    scop_bl_residue_df_domains["pdb_descriptor"] = scop_bl_residue_df_domains["pdb_descriptor"].fillna("")

    scop_bl_domains = scop_bl_residue_df_domains.groupby(scop_cols, dropna = False).agg({"pdb_residue_auth_id": sorted_set, "bound_ligand_auth_id": sorted_set}).reset_index()
    scop_bl_domains["pdb_residue_auth_id"] = scop_bl_domains["pdb_residue_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    scop_bl_domains["bound_ligand_auth_id"] = scop_bl_domains["bound_ligand_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    scop_bl_domains.drop_duplicates(inplace = True)
    scop_bl_domains.to_csv(f"{args.outdir}/scop_bl_domain_ownership.csv", index = False)

    ##ADD BOUND LIGAND TO DFS AND SUGARS

    scop_sugar_residue_df = pd.read_csv(f"{args.scop_sugar_residue_interactions_file}", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    scop_sugar_residue_df.rename(columns = {"ligand_entity_description": "name"}, inplace = True)
    scop_sugar_residue_df["contact_type"] = scop_sugar_residue_df["contact_type"].apply(literal_eval)
    scop_sugar_residue_df_domains = assign_ownership_percentile_categories(scop_sugar_residue_df.copy(), "uniqueID", "scop_id")
    scop_sugar_residue_df_domains["pdb_descriptor"] = scop_sugar_residue_df_domains["pdb_descriptor"].fillna("")
    scop_sugar_domains = scop_sugar_residue_df_domains.groupby(scop_cols, dropna = False).agg({"pdb_residue_auth_id": sorted_set, "bound_ligand_auth_id": sorted_set}).reset_index()
    scop_sugar_domains["pdb_residue_auth_id"] = scop_sugar_domains["pdb_residue_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    scop_sugar_domains["bound_ligand_auth_id"] = scop_sugar_domains["bound_ligand_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    scop_sugar_domains.drop_duplicates(inplace = True)
    scop_sugar_domains.to_csv(f"{args.outdir}/scop_sugar_domain_ownership.csv", index = False)

    scop_combined_domains = pd.concat([scop_bl_domains, scop_sugar_domains])
    scop_combined_domains.to_csv(f"{args.outdir}/scop_combined_domain_ownership.csv", index = False)
    #ADD BOUND LIGAND TO DFS AND SUGARS
    interpro_bl_residue_df = pd.read_csv(f"{args.interpro_bl_residue_interactions_file}", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    interpro_bl_residue_df.rename(columns = {"bound_ligand_name": "name"}, inplace = True)
    interpro_bl_residue_df["contact_type"] = interpro_bl_residue_df["contact_type"].apply(literal_eval)

    interpro_bl_residue_df_d = interpro_bl_residue_df.loc[interpro_bl_residue_df.interpro_type == "Domain"]
    interpro_bl_residue_df_f = interpro_bl_residue_df.loc[interpro_bl_residue_df.interpro_type == "Family"]
    interpro_bl_residue_df_h = interpro_bl_residue_df.loc[interpro_bl_residue_df.interpro_type == "Homologous_superfamily"]

    interpro_bl_residue_df_d_domains = assign_ownership_percentile_categories(interpro_bl_residue_df_d.copy(), "uniqueID", "interpro_accession")
    interpro_bl_residue_df_f_domains = assign_ownership_percentile_categories(interpro_bl_residue_df_f.copy(), "uniqueID", "interpro_accession")
    interpro_bl_residue_df_h_domains = assign_ownership_percentile_categories(interpro_bl_residue_df_h.copy(), "uniqueID", "interpro_accession")

    interpro_bl_domains = pd.concat([interpro_bl_residue_df_d_domains, interpro_bl_residue_df_f_domains, interpro_bl_residue_df_h_domains])

    interpro_cols = ["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords", "protein_entity_id", "protein_entity_ec", "ec_list", "interpro_accession", "interpro_name", "interpro_type", "bm_uniqids", "bound_molecule_display_id", "name", "uniqueID", "bound_ligand_id", "ligand_uniqueID", "type", "total_contact_counts", "domain_contact_counts", "domain_hbond_counts", "domain_contact_perc", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "protein_chain_id", "bound_ligand_struct_asym_id"]
    interpro_bl_domains["pdb_descriptor"] = interpro_bl_domains["pdb_descriptor"].fillna("")
    interpro_bl_domains = interpro_bl_domains.groupby(interpro_cols, dropna = False).agg({"pdb_residue_auth_id": sorted_set, "bound_ligand_auth_id": sorted_set}).reset_index()
    interpro_bl_domains["pdb_residue_auth_id"] = interpro_bl_domains["pdb_residue_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    interpro_bl_domains["bound_ligand_auth_id"] = interpro_bl_domains["bound_ligand_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    interpro_bl_domains.drop_duplicates(inplace = True)
    interpro_bl_domains.to_csv(f"{args.outdir}/interpro_bl_domain_ownership.csv", index = False)


    interpro_sugar_residue_df = pd.read_csv(f"{args.interpro_sugar_residue_interactions_file}", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    interpro_sugar_residue_df.rename(columns = {"ligand_entity_description": "name"}, inplace = True)
    interpro_sugar_residue_df["contact_type"] = interpro_sugar_residue_df["contact_type"].apply(literal_eval)

    interpro_sugar_residue_df_d = interpro_sugar_residue_df.loc[interpro_sugar_residue_df.interpro_type == "Domain"]
    interpro_sugar_residue_df_f = interpro_sugar_residue_df.loc[interpro_sugar_residue_df.interpro_type == "Family"]
    interpro_sugar_residue_df_h = interpro_sugar_residue_df.loc[interpro_sugar_residue_df.interpro_type == "Homologous_superfamily"]

    interpro_sugar_residue_df_d_domains = assign_ownership_percentile_categories(interpro_sugar_residue_df_d.copy(), "uniqueID", "interpro_accession")
    interpro_sugar_residue_df_f_domains = assign_ownership_percentile_categories(interpro_sugar_residue_df_f.copy(), "uniqueID", "interpro_accession")
    interpro_sugar_residue_df_h_domains = assign_ownership_percentile_categories(interpro_sugar_residue_df_h.copy(), "uniqueID", "interpro_accession")

    interpro_sugar_domains = pd.concat([interpro_sugar_residue_df_d_domains,interpro_sugar_residue_df_f_domains,interpro_sugar_residue_df_h_domains])
    interpro_sugar_domains["pdb_descriptor"] = interpro_sugar_domains["pdb_descriptor"].fillna("")
    interpro_sugar_domains = interpro_sugar_domains.groupby(interpro_cols, dropna = False).agg({"pdb_residue_auth_id": sorted_set, "bound_ligand_auth_id": sorted_set}).reset_index()
    interpro_sugar_domains["pdb_residue_auth_id"] = interpro_sugar_domains["pdb_residue_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    interpro_sugar_domains["bound_ligand_auth_id"] = interpro_sugar_domains["bound_ligand_auth_id"].apply(lambda x: "|".join(str(res) for res in x))
    interpro_sugar_domains.drop_duplicates(inplace = True)
    interpro_sugar_domains.to_csv(f"{args.outdir}/interpro_sugar_domain_ownership.csv", index = False)

    interpro_combined_domains = pd.concat([interpro_bl_domains, interpro_sugar_domains])
    interpro_combined_domains.to_csv(f"{args.outdir}/interpro_combined_domain_ownership.csv", index = False)
    
if __name__ == "__main__":
    main()