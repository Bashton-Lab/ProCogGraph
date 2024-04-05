#!/usr/bin/env python
from pdbecif.mmcif_io import CifFileReader
from Bio.ExPASy import Enzyme as EEnzyme
from pathlib import Path
import requests
import os
import numpy as np
from utils import get_terminal_record, get_csdb_from_glycoct, get_glycoct_from_wurcs, get_smiles_from_csdb, process_ec_records, Neo4jConnection
import json
from urllib.parse import quote
import pandas as pd
import xml.etree.ElementTree as ET
import re
import argparse
from bs4 import BeautifulSoup
import yaml
import concurrent.futures
from rich.progress import Progress
from rich.console import Console
import pickle
import gzip 
import re

def return_partial_EC_list(ec, total_ec_list):
    if not isinstance(ec, str) and np.isnan(ec):
        return np.nan
    elif "-" in ec:
        replacement_character = r'.'
        modified_ec = re.sub(r'\.', r"_", ec)
        modified_ec = modified_ec.replace("-", ".")
        total_ec_list = [re.sub(r'\.', r"_", item) for item in total_ec_list]
        # Use re.match() to check if the modified string matches any item in the match_list
        matching_ec = [ec for ec in total_ec_list if re.match(modified_ec, ec)]
        matching_ec = [re.sub(r'_', r".", item) for item in matching_ec]
        return(matching_ec)
    else:
        return [ec]

def get_updated_enzyme_records(df, ec_records_df, ec_col = "protein_entity_ec"):
    ec_list = ec_records_df.ID.unique() ##fill the partial ec records using the original ec ids from the expasy enzyme list
    
    residue_ec_records = df[[ec_col]].drop_duplicates()
    residue_ec_records["protein_entity_ec_copy"] = residue_ec_records[ec_col]
    residue_ec_records["protein_entity_ec_copy"] = residue_ec_records.protein_entity_ec_copy.str.split(",")
    residue_ec_records = residue_ec_records.explode("protein_entity_ec_copy")
    residue_ec_records["protein_entity_ec_copy"] = residue_ec_records.protein_entity_ec_copy.str.strip()
    residue_ec_records["ec_list"] = residue_ec_records.protein_entity_ec_copy.apply(lambda x: return_partial_EC_list(x, ec_list))
    residue_ec_records = residue_ec_records.explode("ec_list")
    residue_ec_records = residue_ec_records.merge(ec_records_df[["ID", "TRANSFER"]], left_on = "ec_list", right_on = "ID", how = "left")
    residue_ec_records["TRANSFER"] = residue_ec_records["TRANSFER"].fillna("")

    # anythin with NAN now in ID/transfer doesnt actually exist in the expasy enzyme list - so is incorrect.

    residue_ec_records_grouped = residue_ec_records.groupby(ec_col).agg({"TRANSFER": set}).reset_index()
    residue_ec_records_grouped["TRANSFER"] = residue_ec_records_grouped["TRANSFER"].apply(lambda x: ",".join(x) if x != "" else "")
    residue_ec_records_grouped.rename(columns = {"TRANSFER" : "ec_list"}, inplace = True)
    
    df_merged = df.merge(residue_ec_records_grouped, on = ec_col, how = "left", indicator = True)
    assert(len(df_merged.loc[df_merged["_merge"] != "both"]) == 0)
    df_merged.drop(columns = "_merge", inplace = True)
    df_merged = df_merged.loc[df_merged.ec_list != ""] #remove any rows where the ec_list is empty - we cant process these anyway.
    return(df_merged)

def clean_and_merge_scop_col(df, column_id, description_df):
    level = df[column_id].str.split("=").str.get(0).values[0]
    df[column_id] = df[column_id].str.split("=").str.get(1).astype(int)
    df = df.merge(description_df.loc[description_df.level == level, ["level_sunid", "level", "level_description"]],left_on = column_id, right_on = "level_sunid", indicator = True)
    df.rename(columns = {"level_description": f"{level}_description"}, inplace = True)
    assert len(df.loc[df._merge != "both"]) == 0
    df.drop(columns = ["_merge", "level_sunid", "level"], inplace = True)
    return df

def process_row(conn, progress, task, pdb, query):
    result = conn.query(query, db='neo4j', pdb_id=pdb)
    progress.update(task, advance=1)
    return result

def process_row_sifts(conn, progress, task, row, query):
    result = conn.query(query, db='neo4j', pdb_id=row.PDB, chain=row.CHAIN)
    progress.update(task, advance=1)
    return result

def parse_table_data(elem):
    data = {}
    for row in elem.findall('row'):
        row_data = {}
        for field in row.findall('field'):
            row_data[field.attrib['name']] = field.text
        data[len(data)] = row_data
    return data

def assign_ownership_percentile_categories(ligands_df, unique_id = "uniqueID", domain_grouping_key = "cath_domain"):
    ligands_df["total_contact_counts"]  =  ligands_df.groupby([unique_id])["contact_type_count"].transform("sum")
    ligands_df[["domain_contact_counts", "domain_hbond_counts", "domain_covalent_counts"]]  = ligands_df.groupby([unique_id, domain_grouping_key])[["contact_type_count", "hbond_count", "covalent_count"]].transform("sum")
    ligands_df["domain_hbond_perc"] = ligands_df.domain_hbond_counts / ligands_df.total_contact_counts
    ligands_df["domain_contact_perc"] = ligands_df.domain_contact_counts / ligands_df.total_contact_counts
    ligands_df["num_non_minor_domains"] = ligands_df.groupby([unique_id])["domain_contact_perc"].transform(lambda x: len(x[x > 0.1]))
    ligands_df["domain_ownership"] = np.where(
        ligands_df["domain_contact_perc"] == 1, "unique",
        np.where(
            ligands_df["domain_contact_perc"] >= 0.9, "dominant",
            np.where(
                (ligands_df["domain_contact_perc"] >= 0.5)
                & (ligands_df["domain_contact_perc"] < 0.9) & (ligands_df["num_non_minor_domains"] == 1), "major",
                np.where(
                (ligands_df["domain_contact_perc"] >= 0.5)
                & (ligands_df["domain_contact_perc"] < 0.9) & (ligands_df["num_non_minor_domains"] >1), "major_partner",
                    np.where(
                    (ligands_df["domain_contact_perc"] >= 0.1)
                    & (ligands_df["domain_contact_perc"] < 0.5) & (ligands_df["num_non_minor_domains"] >1), "partner",
                        np.where(
                        ligands_df["domain_contact_perc"] < 0.1, "minor", np.nan)
                    )
                )
            )
        )
    )
    
    return ligands_df

def extract_domain_annotations(xml_file):
    with gzip.open(xml_file, 'rb') as f:
        tree = ET.parse(f)
        root = tree.getroot()

        # Find all interpro elements in the XML file
        interpro_elements = root.findall('.//interpro')

        pfam_annotations = {}
        superfamily_annotations = {}
        gene3d_annotations = {}
        # Iterate through each interpro element
        for interpro in interpro_elements:
            interpro_id = interpro.attrib['id']
            # Find db_xref elements with db attribute as PFAM
            pfam_refs = interpro.findall('.//db_xref[@db="PFAM"]')
            superfamily_refs = interpro.findall('.//db_xref[@db="SSF"]')
            gene3d_refs = interpro.findall('.//db_xref[@db="CATHGENE3D"]')
            pfam_accessions = []
            superfamily_accessions = []
            gene3d_accessions = []
            # Extract PFAM accessions for the current interpro element
            for pfam_ref in pfam_refs:
                pfam_accessions.append("PFAM:" + pfam_ref.attrib.get('dbkey'))
            for superfamily_ref in superfamily_refs:
                superfamily_accessions.append("SUPERFAMILY:" + superfamily_ref.attrib.get('dbkey'))
            for gene3d_ref in gene3d_refs:
                gene3d_accessions.append(gene3d_ref.attrib.get('dbkey')) #no prefix for gene3d as it is prefixed in ref

            # Store PFAM annotations for the interpro ID
            pfam_annotations[interpro_id] = pfam_accessions
            superfamily_annotations[interpro_id] = superfamily_accessions
            gene3d_annotations[interpro_id] = gene3d_accessions

    interpro_annotations = pd.DataFrame([pfam_annotations, superfamily_annotations, gene3d_annotations], index = ["pfam_annotations", "superfamily_annotations", "gene3d_annotations"]).T
    interpro_annotations["dbxref"] = interpro_annotations["pfam_annotations"].str.join("|") + "|" + interpro_annotations["superfamily_annotations"].str.join("|") + "|" + interpro_annotations["gene3d_annotations"].str.join("|")
    interpro_annotations["dbxref"] = interpro_annotations["dbxref"].str.rstrip("|").str.lstrip("|")
    interpro_annotations["dbxref"] = interpro_annotations["dbxref"].replace("", np.nan)
    return interpro_annotations[["dbxref"]]

def main():

    parser = argparse.ArgumentParser(description = 'TO DO')
    parser.add_argument('--neo4j_bolt_uri', default = 'bolt://localhost:7687', type = str,
        help = "")
    parser.add_argument('--neo4j_user', default = 'neo4j', type = str,
        help = "")
    parser.add_argument('--neo4j_password', type = str,
        help = "")
    parser.add_argument('--outdir', type = str, default = "parity_calcs",
        help = ""),
    parser.add_argument('--enzyme_dat_file', type = str, default = "enzyme.dat",
        help = ""),
    parser.add_argument('--enzyme_class_file', type = str, default = "enzclass.txt",
        help = "")
    parser.add_argument('--pdbe_graph_yaml', type = str, default = "pdbe_graph.yaml",
        help = "")
    parser.add_argument('--glycoct_cache', type = str,
        help = "glycoct cache file")
    parser.add_argument('--smiles_cache', type = str,
        help = "smiles cache file")
    parser.add_argument('--csdb_linear_cache', type = str,
        help = "csdb linear cache file")
    parser.add_argument('--sifts_ec_mapping', type = str,
        help = "sifts ec mapping file")
    parser.add_argument('--threads', type = int, default = 1,
        help = "Number of threads to use")
    parser.add_argument('--scop_domains_info_file', type = str,
        help = "scop domains info file")
    parser.add_argument('--scop_descriptions_file', type = str,
        help="scop descriptions file")
    parser.add_argument('--interpro_xml', metavar='interpro_xml', type=str,
        help = "path to interpro xml file (gzipped)")
    parser.add_argument('--pfam_clan_rels', type = str,
        help = "pfam clan relationships file clan_membership.txt.gz")
    parser.add_argument('--pfam_clans', type = str,
        help = "pfam clans file clan.txt.gz")
    parser.add_argument('--pdbecif_dir', type = str, default = "pdbecif",
        help = "path to directory containing cif files")
    parser.add_argument('--domain_contact_cutoff', type = int, default = 3,
        help = "minimum number of contacts for a domain to be considered in the analysis (default 3)")
    
    
    args = parser.parse_args()

    console = Console()
    Path(f"{args.outdir}").mkdir(parents=True, exist_ok=True)

    if args.glycoct_cache:
        glycoct_cache = pd.read_pickle(args.glycoct_cache)
    else:
        glycoct_cache = pd.DataFrame(columns = ["WURCS", "glycoct"])
    if args.smiles_cache:
        smiles_cache = pd.read_pickle(args.smiles_cache)
    else:
        smiles_cache = pd.DataFrame(columns = ["csdb", "descriptor"])
    if args.csdb_linear_cache:
        csdb_linear_cache = pd.read_pickle(args.csdb_linear_cache)
    else:
        csdb_linear_cache = pd.DataFrame(columns = ["glycoct", "csdb"])
        
    #from https://stackoverflow.com/questions/1773805/how-can-i-parse-a-yaml-file-in-python
    with open(f"{args.pdbe_graph_yaml}", "r") as stream:
        try:
            pdbe_graph_queries = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    bl_queries = {
        "CATH" : {"query" : pdbe_graph_queries["cath_bl_query"], "domain_id": "cath_domain"},
        "SCOP" : {"query" : pdbe_graph_queries["scop_bl_query"], "domain_id": "scop_id"},
        "PFAM" : {"query" : pdbe_graph_queries["pfam_bl_query"], "domain_id": "pfam_accession"},
        "InterProHomologousSuperfamily" : {"query": pdbe_graph_queries["interpro_h_bl_query"], "domain_id": "interpro_accession"}}

    bs_queries = {
        "CATH" : {"query" : pdbe_graph_queries["cath_sugar_query"], "domain_id": "cath_domain"},
        "SCOP" : {"query" : pdbe_graph_queries["scop_sugar_query"], "domain_id": "scop_id"},
        "PFAM" : {"query" : pdbe_graph_queries["pfam_sugar_query"], "domain_id": "pfam_accession"},
        "InterProHomologousSuperfamily" : {"query" : pdbe_graph_queries["interpro_h_sugar_query"], "domain_id": "interpro_accession"}}

    console.print("Connecting to neo4j")
    conn = Neo4jConnection(uri=args.neo4j_bolt_uri, user=args.neo4j_user, pwd=args.neo4j_password)
    console.print("Connected to neo4j")
    console.print("Generating EC record dataframe")
    
    with Progress() as progress:
            
        ec_records_df_grouped = process_ec_records(args.enzyme_dat_file , args.enzyme_class_file)
        ec_records_df = ec_records_df_grouped.explode("ID")
        sifts_chains = pd.read_csv(f"{args.sifts_ec_mapping}", sep = "\t", comment="#")
        sifts_chains_ec = sifts_chains.loc[sifts_chains.EC_NUMBER != "?"]
        sifts_chains_ec = sifts_chains_ec.groupby(["PDB", "CHAIN"]).agg({"ACCESSION": set, "EC_NUMBER": set}).reset_index() #some chains have multiple uniprot accessions, we group these into a list along with their associated ec's
        sifts_chains_ec["EC_NUMBER"] = sifts_chains_ec["EC_NUMBER"].apply(lambda x: ",".join(x)) #join the list of EC numbers into a single string for ec records function 
        sifts_chains_ec["ACCESSION"] = sifts_chains_ec["ACCESSION"].apply(lambda x: "|".join(x)) #join the list of uniprot accessions with a pipe for downstream neo4j integration
        sifts_chains_ec = get_updated_enzyme_records(sifts_chains_ec, ec_records_df, ec_col = "EC_NUMBER")
        sifts_chains_ec.rename(columns = {"EC_NUMBER": "protein_entity_ec", "ACCESSION" : "uniprot_accession"}, inplace = True)
        

        all_pdbs = sifts_chains_ec.PDB.unique()
        total_rows = len(all_pdbs)

        bl_results = {}
        bl_results_unmatched = {}
        bs_results = {}
        bs_results_unmatched = {}

        scop_domains_info = pd.read_csv(f"{args.scop_domains_info_file}", sep = "\t", comment = "#", header = None, names = ["scop_id", "pdb_id", "scop_description", "sccs", "domain_sunid", "ancestor_sunid"])
        scop_id_levels = ["cl_id", "cf_id", "sf_id", "fa_id", "dm_id", "sp_id", "px_id"]
        scop_domains_info[scop_id_levels] = scop_domains_info.ancestor_sunid.str.split(",", expand = True)
        scop_descriptions = pd.read_csv(f"{args.scop_descriptions_file}", sep = "\t", comment = "#" , header = None, names = ["level_sunid", "level", "level_sccs", "level_sid", "level_description"])

        for column in scop_id_levels:
            scop_domains_info = clean_and_merge_scop_col(scop_domains_info, column, scop_descriptions)
        
        scop_domains_info.drop(columns = ["pdb_id", "scop_description"], inplace = True)

        class_codes = scop_domains_info[["cl_id", "cl_description"]].drop_duplicates()
        fold_codes = scop_domains_info[["cf_id", "cf_description"]].drop_duplicates()
        superfamily_codes = scop_domains_info[["sf_id", "sf_description"]].drop_duplicates()

        interpro_annotations = extract_domain_annotations(args.interpro_xml)

        pfam_clan_rels = pd.read_csv(f"{args.pfam_clan_rels}", sep = "\t", header = None, names = ["clan", "pfam"])
        pfam_clans = pd.read_csv(f"{args.pfam_clans}", sep = "\t", comment = "#", header = None, names = ["clan_acc", "clan_id", "previous_id", "clan_description", "clan_author", "deposited_by", "clan_comment", "updated", "created", "version", "number_structures", "number_archs", "number_species", "number_sequences", "competed", "uniprot_competed"])
        pfam_df = pfam_clan_rels.merge(pfam_clans[["clan_acc", "clan_description", "clan_comment"]], left_on = "clan", right_on = "clan_acc", how = "left", indicator = True)
        assert(len(pfam_df.loc[pfam_df._merge != "both"]) == 0)
        pfam_df.drop(columns = "_merge", inplace = True)

        if not os.path.exists(f"{args.outdir}/bl_results.pkl"):

            for db, data in bl_queries.items():
                query = data["query"]
                domain_identifier = data["domain_id"]
                task = progress.add_task(f"[cyan]Processing {db} bound ligands...", total=total_rows)
                results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
                    futures = [executor.submit(process_row, conn = conn, progress = progress, task = task, pdb = pdb, query = query) for pdb in all_pdbs]
                    for future in concurrent.futures.as_completed(futures):
                        if future.result() != None:
                            results.extend(future.result())

                # Convert results to DataFrame
                result_df = pd.DataFrame([dict(_) for _ in results])
                result_df_contact_filtered = result_df.loc[result_df.contact_type_count >= args.domain_contact_cutoff].copy()
                result_df_ec = result_df_contact_filtered.merge(sifts_chains_ec, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left", indicator = True) #keeping only pdbs with sifts ec annotations
                result_df_ec_unmatched = result_df_ec.loc[result_df_ec._merge != "both"].copy().drop(columns = ["_merge", "PDB", "CHAIN"])
                result_df_ec = result_df_ec.loc[result_df_ec._merge == "both"].copy().drop(columns = ["_merge", "PDB", "CHAIN"])
                
                if db == "SCOP":
                    result_df_ec = result_df_ec.merge(scop_domains_info, how = "left", on = "scop_id", indicator = True)
                    assert(len(result_df_ec.loc[result_df_ec._merge != "both"]) == 0)
                    result_df_ec.drop(columns = "_merge", inplace = True)
                elif db == "PFAM":
                    result_df_ec = result_df_ec.merge(pfam_df, left_on = "pfam_accession", right_on = "pfam", how = "left")
                elif db == "InterProHomologousSuperfamily":
                    result_df_ec = result_df_ec.merge(interpro_annotations, left_on = "interpro_accession", right_index = True, how = "left")
                    result_df_ec = result_df_ec.loc[(result_df_ec.dbxref.isna() == False) & ((result_df_ec.dbxref.str.contains("SUPERFAMILY")) | (result_df_ec.dbxref.str.contains("G3DSA")))] 

                console.print("Assigning ownership categories")
                result_df_ec_ownership = assign_ownership_percentile_categories(result_df_ec, unique_id = "uniqueID", domain_grouping_key = domain_identifier)
                bl_results[db] = result_df_ec_ownership
                bl_results_unmatched[db] = result_df_ec_unmatched
            del(results)
            console.print("Saving bound ligand results to file")
            with open(f"{args.outdir}/bl_results.pkl", 'wb') as f:
                pickle.dump(bl_results, f)
            with open(f"{args.outdir}/bl_results_sifts_unmatched.pkl", 'wb') as f:
                pickle.dump(bl_results_unmatched, f)
        else:
            with open(f"{args.outdir}/bl_results.pkl", 'rb') as f:
                bl_results = pickle.load(f)
            console.print(f"Loaded bound ligand results from file {args.outdir}/bl_results.pkl")

        if not os.path.exists(f"{args.outdir}/bs_results.pkl"):
            for db, data in bs_queries.items():
                query = data["query"]
                domain_identifier = data["domain_id"]
                task = progress.add_task(f"[cyan]Processing {db} bound sugars...", total=total_rows)
                results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
                    futures = [executor.submit(process_row, conn = conn, progress = progress, task = task, pdb = pdb, query = query) for pdb in all_pdbs]
                    for future in concurrent.futures.as_completed(futures):
                        if future.result() != None:
                            results.extend(future.result())

                # Convert results to DataFrame
                result_df = pd.DataFrame([dict(_) for _ in results])
                result_df_contact_filtered = result_df.loc[result_df.contact_type_count >= args.domain_contact_cutoff].copy()
                result_df["ligand_entity_id_numerical"] = result_df["ligand_entity_id_numerical"].astype(int)
                result_df_ec = result_df.merge(sifts_chains_ec, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left", indicator = True) #keeping only pdbs with sifts ec annotations
                result_df_ec_unmatched = result_df_ec.loc[result_df_ec._merge != "both"].copy().drop(columns = ["_merge", "PDB", "CHAIN"])
                result_df_ec = result_df_ec.loc[result_df_ec._merge == "both"].copy().drop(columns = ["_merge", "PDB", "CHAIN"]) 
                if db == "SCOP":
                    result_df_ec = result_df_ec.merge(scop_domains_info, how = "left", on = "scop_id", indicator = True)
                    assert(len(result_df_ec.loc[result_df_ec._merge != "both"]) == 0)
                    result_df_ec.drop(columns = "_merge", inplace = True)
                elif db == "PFAM":
                    result_df_ec = result_df_ec.merge(pfam_df, left_on = "pfam_accession", right_on = "pfam", how = "left")
                elif db == "InterProHomologousSuperfamily":
                    result_df_ec = result_df_ec.merge(interpro_annotations, left_on = "interpro_accession", right_index = True, how = "left")
                    result_df_ec = result_df_ec.loc[(result_df_ec.dbxref.isna() == False) & ((result_df_ec.dbxref.str.contains("SUPERFAMILY")) | (result_df_ec.dbxref.str.contains("G3DSA")))] 

                console.print("Assigning ownership categories")
                result_df_ec_ownership = assign_ownership_percentile_categories(result_df_ec, unique_id = "uniqueID", domain_grouping_key = domain_identifier)
                bs_results[db] = result_df_ec
                bs_results_unmatched[db] = result_df_ec_unmatched
            del(results)
            console.print("Saving bound sugar results to file")
            with open(f"{args.outdir}/bs_results.pkl", 'wb') as f:
                pickle.dump(bs_results, f)
            with open(f"{args.outdir}/bs_results_sifts_unmatched.pkl", 'wb') as f:
                pickle.dump(bs_results_unmatched, f)
        else:
            with open(f"{args.outdir}/bs_results.pkl", 'rb') as f:
                bs_results = pickle.load(f)
            console.print(f"Loaded bound sugar results from file {args.outdir}/bs_results.pkl")

    cath_pdb_residue_interactions_bl = bl_results["CATH"]
    scop_pdb_residue_interactions_bl = bl_results["SCOP"]
    pfam_pdb_residue_interactions_bl = bl_results["PFAM"]
    interpro_pdb_residue_interactions_bl = bl_results["InterProHomologousSuperfamily"]
    
    #check what we are using this for again
    bound_molecules_ligands = pd.concat([cath_pdb_residue_interactions_bl[["bm_ids", "bound_ligand_id"]].drop_duplicates(), scop_pdb_residue_interactions_bl[["bm_ids", "bound_ligand_id"]].drop_duplicates(),
            interpro_pdb_residue_interactions_bl[["bm_ids", "bound_ligand_id"]].drop_duplicates()])
    
    bound_molecules_ligands.to_csv(f"{args.outdir}/bound_molecules_ligands.csv.gz", index = False, compression = "gzip")
    
    cath_pdb_residue_interactions_bs = bs_results["CATH"]
    scop_pdb_residue_interactions_bs = bs_results["SCOP"]
    pfam_pdb_residue_interactions_bs = bs_results["PFAM"]
    interpro_pdb_residue_interactions_bs = bs_results["InterProHomologousSuperfamily"]
     
    if not os.path.exists(f"{args.outdir}/bound_molecules_sugars_wurcs.csv.gz"):
        bound_molecules_sugars = pd.concat([
                cath_pdb_residue_interactions_bs[["pdb_id", "bm_ids", "ligand_entity_id", "uniqueID", "description", "ligand_entity_id_numerical", "protein_entity_ec", "ec_list"]].drop_duplicates() , 
                scop_pdb_residue_interactions_bs[["pdb_id", "bm_ids", "ligand_entity_id", "uniqueID", "description", "ligand_entity_id_numerical", "protein_entity_ec", "ec_list"]].drop_duplicates() ,
                pfam_pdb_residue_interactions_bs[["pdb_id", "bm_ids", "ligand_entity_id", "uniqueID", "description", "ligand_entity_id_numerical", "protein_entity_ec", "ec_list"]].drop_duplicates() ,
                interpro_pdb_residue_interactions_bs[["pdb_id", "bm_ids", "ligand_entity_id", "uniqueID", "description", "ligand_entity_id_numerical", "protein_entity_ec", "ec_list"]].drop_duplicates()
                ])

        cif_ids = bound_molecules_sugars.pdb_id.unique()
        cif_dir = f"{args.pdbecif_dir}"
        Path(f"{cif_dir}").mkdir(parents=True, exist_ok=True)
        task = progress.add_task(f"[cyan]Downloading bound sugars CIF files...", total=len(cif_ids))
        
        for cif_id in cif_ids:
            cif_path = f"{cif_dir}/{cif_id}_updated.cif"
            if not os.path.exists(cif_path):
                response = requests.get(f'http://www.ebi.ac.uk/pdbe/entry-files/download/{cif_id}_updated.cif')

                with open(cif_path, 'wb') as fp:
                    fp.write(response.content)
            progress.update(task, advance=1)


        branched_entity_list = []
        reader = CifFileReader()
        task = progress.add_task(f"[cyan]Extracting WURCS from bound sugars CIF files...", total=len(cif_ids))
        for cif_id in cif_ids:
            cif_path = f'{cif_dir}/{cif_id}_updated.cif'
            cif_dict = reader.read(cif_path, output='cif_dictionary')
            cif_df = pd.DataFrame(cif_dict)
            if "_pdbx_entity_branch_descriptor" in cif_df.index:
                branch_data = cif_df.loc["_pdbx_entity_branch_descriptor"].values[0]
                if type(next(iter(branch_data.values()))) == str:  # If only one descriptor
                    branched_entity_descriptors = pd.DataFrame({key: [value] for key, value in branch_data.items()})  # Convert scalar value to DataFrame
                else:  # If multiple descriptors
                    branched_entity_descriptors = pd.DataFrame(branch_data)
                branched_entity_descriptors["pdb_id"] = cif_id

                branched_entity_list.append(branched_entity_descriptors)
            progress.update(task, advance=1)
        branched_entity_df = pd.concat(branched_entity_list)
        branched_entity_df["entity_id"] = branched_entity_df.entity_id.astype("int")
        branched_entity_df.loc[branched_entity_df.type == "WURCS"].pdb_id.nunique()
        
        sugar_wurcs = branched_entity_df.loc[branched_entity_df.type == "WURCS"].groupby(["pdb_id", "entity_id"]).head(1).reset_index()

        sugar_wurcs.rename(columns = {"descriptor" : "WURCS"}, inplace = True)
        bound_molecules_sugars_wurcs = bound_molecules_sugars.merge(sugar_wurcs[["entity_id", "WURCS", "pdb_id"]], left_on = ["pdb_id", "ligand_entity_id_numerical"], right_on = ["pdb_id", "entity_id"], how = "left", indicator = True)
        
        bound_molecules_sugars_wurcs.loc[bound_molecules_sugars_wurcs._merge != "both", "WURCS"] = "WURCS not available"
        bound_molecules_sugars_wurcs.drop(columns = ["_merge"], inplace = True)

        
        bound_molecules_sugars_wurcs["ec_list"] = bound_molecules_sugars_wurcs.ec_list.str.split(",")
        bound_molecules_sugars_wurcs = bound_molecules_sugars_wurcs.explode("ec_list")
        bound_molecules_sugars_wurcs.drop(columns = "protein_entity_ec", inplace = True)
        bound_molecules_sugars_wurcs.to_csv(f"{args.outdir}/bound_molecules_sugars_wurcs.csv.gz", compression = "gzip", index = False)
    else:
        bound_molecules_sugars_wurcs = pd.read_csv(f"{args.outdir}/bound_molecules_sugars_wurcs.csv.gz", compression = "gzip")
        console.print(f"Loaded bound sugar WURCS from file {args.outdir}/bound_molecules_sugars_wurcs.csv.gz")
    
    if not os.path.exists(f"{args.outdir}/bound_sugars_to_score.pkl") or not os.path.exists(f"{args.outdir}/bound_molecules_sugars_smiles.pkl") or not os.path.exists(f"{args.outdir}/cath_pdb_residue_interactions_bs_index.pkl"):
        bound_sugars_to_score = bound_molecules_sugars_wurcs.loc[bound_molecules_sugars_wurcs.WURCS != "WURCS not available", ["description","WURCS", "ec_list"]].drop_duplicates()
        bound_sugars_to_score = bound_sugars_to_score.groupby(["description","WURCS"]).agg({"ec_list": set}).reset_index()

        bound_sugars_to_score["glycoct"] = bound_sugars_to_score["WURCS"].apply(lambda x: get_glycoct_from_wurcs(x, glycoct_cache))
        new_glycoct_values = bound_sugars_to_score.loc[bound_sugars_to_score.WURCS.isin(glycoct_cache.WURCS.values) == False, ["glycoct","WURCS"]].drop_duplicates()
        glycoct_cache = pd.concat([glycoct_cache, new_glycoct_values], ignore_index = True)
        glycoct_cache.to_pickle(f"{args.glycoct_cache}")
        bound_sugars_to_score = bound_sugars_to_score.loc[bound_sugars_to_score.glycoct.isna() == False]

        bound_sugars_to_score["csdb"] = bound_sugars_to_score["glycoct"].apply(lambda x: get_csdb_from_glycoct(x, csdb_linear_cache))
        new_csdb_values = bound_sugars_to_score.loc[bound_sugars_to_score.glycoct.isin(csdb_linear_cache.glycoct.values) == False, ["csdb","glycoct"]].drop_duplicates()
        csdb_linear_cache = pd.concat([csdb_linear_cache, new_csdb_values], ignore_index = True)
        csdb_linear_cache.to_pickle(f"{args.csdb_linear_cache}")
        bound_sugars_to_score = bound_sugars_to_score.loc[bound_sugars_to_score.csdb.isna() == False]
        bound_sugars_to_score["descriptor"] = bound_sugars_to_score["csdb"].apply(lambda x: get_smiles_from_csdb(x, smiles_cache))
        new_smiles_values = bound_sugars_to_score.loc[bound_sugars_to_score.csdb.isin(smiles_cache.csdb.values) == False, ["descriptor","csdb"]].drop_duplicates()
        smiles_cache = pd.concat([smiles_cache, new_smiles_values], ignore_index = True)
        smiles_cache.to_pickle(f"{args.smiles_cache}")
        bound_sugars_to_score = bound_sugars_to_score.loc[bound_sugars_to_score.descriptor.isna() == False]

        bound_sugars_to_score = bound_sugars_to_score.reset_index()
        bound_sugars_to_score.drop(columns = ["index"], inplace = True)
        bound_sugars_to_score = bound_sugars_to_score.reset_index().rename(columns = {"index": "ligand_index"})

        bound_molecules_sugars_smiles = bound_molecules_sugars_wurcs.merge(bound_sugars_to_score[["description", "ligand_index", "WURCS", "descriptor"]], on = ["description","WURCS"], how = "left")

        bound_sugars_to_score["bl_name"] = bound_sugars_to_score["description"]

        bound_sugars_to_score.rename(columns = {"ligand_index": "ligand_entity_id"}, inplace = True) #do this to run sugars in parity calcs
        bound_sugars_to_score.to_pickle(f"{args.outdir}/bound_sugars_to_score.pkl")
        
        missing_ligand_index = bound_molecules_sugars_smiles.loc[bound_molecules_sugars_smiles.descriptor.isna(), ["uniqueID"]].drop_duplicates()
        missing_ligand_index["missing_ligand_index"] = missing_ligand_index.reset_index(drop=True).reset_index().index + bound_molecules_sugars_smiles.ligand_index.max() + 1

        bound_molecules_sugars_smiles = bound_molecules_sugars_smiles.merge(missing_ligand_index, on = ["uniqueID"], how = "left")
        bound_molecules_sugars_smiles["ligand_index"].fillna(bound_molecules_sugars_smiles["missing_ligand_index"], inplace=True)
        bound_molecules_sugars_smiles.drop(columns = "missing_ligand_index", inplace = True)
        bound_molecules_sugars_smiles["descriptor"].fillna("SMILES unavailable", inplace = True)
        bound_molecules_sugars_smiles.to_pickle(f"{args.outdir}/bound_molecules_sugars_smiles.pkl")

        cath_pre_shape = len(cath_pdb_residue_interactions_bs)
        cath_pdb_residue_interactions_bs_index = cath_pdb_residue_interactions_bs.merge(bound_molecules_sugars_smiles[["ligand_index", "uniqueID", "descriptor"]].drop_duplicates().rename(columns = {"ligand_index": "ligand_uniqueID"}), on = "uniqueID", how = "left", indicator = True, validate = "many_to_one")
        assert(len(cath_pdb_residue_interactions_bs_index) == cath_pre_shape)
        assert(len(cath_pdb_residue_interactions_bs_index.loc[cath_pdb_residue_interactions_bs_index._merge != "both"]) == 0)
        cath_pdb_residue_interactions_bs_index.drop(columns = ["_merge"], inplace = True)
        cath_pdb_residue_interactions_bs_index.to_csv(f"{args.outdir}/cath_pdb_residue_interactions_bs_index.csv.gz", index = False, compression = "gzip")

        scop_pre_shape = len(scop_pdb_residue_interactions_bs)
        scop_pdb_residue_interactions_bs_index = scop_pdb_residue_interactions_bs.merge(bound_molecules_sugars_smiles[["ligand_index", "uniqueID", "descriptor"]].drop_duplicates().rename(columns = {"ligand_index": "ligand_uniqueID"}), on = "uniqueID", how = "left", indicator = True, validate = "many_to_one")
        assert(len(scop_pdb_residue_interactions_bs_index) == scop_pre_shape)
        assert(len(scop_pdb_residue_interactions_bs_index.loc[scop_pdb_residue_interactions_bs_index._merge != "both"]) == 0)
        scop_pdb_residue_interactions_bs_index.drop(columns = ["_merge"], inplace = True)
        scop_pdb_residue_interactions_bs_index.to_csv(f"{args.outdir}/scop_pdb_residue_interactions_bs_index.csv.gz", index = False, compression = "gzip")

        pfam_pre_shape = len(pfam_pdb_residue_interactions_bs)
        pfam_pdb_residue_interactions_bs_index = pfam_pdb_residue_interactions_bs.merge(bound_molecules_sugars_smiles[["ligand_index", "uniqueID", "descriptor"]].drop_duplicates().rename(columns = {"ligand_index": "ligand_uniqueID"}), on = "uniqueID", how = "left", indicator = True, validate = "many_to_one")
        assert(len(pfam_pdb_residue_interactions_bs_index) == pfam_pre_shape)
        assert(len(pfam_pdb_residue_interactions_bs_index.loc[pfam_pdb_residue_interactions_bs_index._merge != "both"]) == 0)
        pfam_pdb_residue_interactions_bs_index.drop(columns = ["_merge"], inplace = True)
        pfam_pdb_residue_interactions_bs_index.to_csv(f"{args.outdir}/pfam_pdb_residue_interactions_bs_index.csv.gz", index = False, compression = "gzip")

        interpro_pre_shape = len(interpro_pdb_residue_interactions_bs)
        interpro_pdb_residue_interactions_bs_index = interpro_pdb_residue_interactions_bs.merge(bound_molecules_sugars_smiles[["ligand_index", "uniqueID", "descriptor"]].drop_duplicates().rename(columns = {"ligand_index": "ligand_uniqueID"}), on = "uniqueID", how = "left", indicator = True, validate = "many_to_one")
        assert(len(interpro_pdb_residue_interactions_bs_index) == interpro_pre_shape)
        assert(len(interpro_pdb_residue_interactions_bs_index.loc[interpro_pdb_residue_interactions_bs_index._merge != "both"]) == 0)
        interpro_pdb_residue_interactions_bs_index.drop(columns = ["_merge"], inplace = True)
        interpro_pdb_residue_interactions_bs_index.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_bs_index.csv.gz", index = False, compression = "gzip")
        
    else:
        bound_molecules_sugars_smiles = pd.read_pickle(f"{args.outdir}/bound_molecules_sugars_smiles.pkl")
        bound_sugars_to_score = pd.read_pickle(f"{args.outdir}/bound_sugars_to_score.pkl")
        cath_pdb_residue_interactions_bs_index = pd.read_csv(f"{args.outdir}/cath_pdb_residue_interactions_bs_index.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        scop_pdb_residue_interactions_bs_index = pd.read_csv(f"{args.outdir}/scop_pdb_residue_interactions_bs_index.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        pfam_pdb_residue_interactions_bs_index = pd.read_csv(f"{args.outdir}/pfam_pdb_residue_interactions_bs_index.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        interpro_pdb_residue_interactions_bs_index = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_bs_index.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

        console.print(f"Loaded bound molecules sugars smiles from file {args.outdir}/bound_molecules_sugars_smiles.pkl and bound sugars to score from {args.outdir}/bound_sugars_to_score.pkl and bound sugar residue interaction ID files")
    
    
    if not os.path.exists(f"{args.outdir}/bound_ligands_to_score.pkl"):
        bound_ligands_to_score = pd.concat([cath_pdb_residue_interactions_bl[["description", "bound_ligand_name", "descriptor", "ec_list"]],
                                            scop_pdb_residue_interactions_bl[["description", "bound_ligand_name", "descriptor", "ec_list"]],
                                            pfam_pdb_residue_interactions_bl[["description", "bound_ligand_name", "descriptor", "ec_list"]],
                                            interpro_pdb_residue_interactions_bl[["description", "bound_ligand_name", "descriptor", "ec_list"]]]).drop_duplicates()

        bound_ligands_to_score = bound_ligands_to_score.groupby(["bound_ligand_name", "descriptor"]).agg({"ec_list": set, "description": "first"}).reset_index()

        bound_ligands_to_score = bound_ligands_to_score.reset_index().rename(columns = {"index" : "ligand_entity_id"})
        bound_ligands_to_score["ligand_entity_id"] = bound_ligands_to_score["ligand_entity_id"] + bound_molecules_sugars_smiles.ligand_index.max() + 1 #plus one because of 0 index to avoid overlaps

        cath_pdb_residue_interactions_bl_id = cath_pdb_residue_interactions_bl.merge(bound_ligands_to_score[["bound_ligand_name", "ligand_entity_id"]].rename(columns = {"ligand_entity_id": "ligand_uniqueID"}), on = "bound_ligand_name", how = "left", indicator = True)
        assert(len(cath_pdb_residue_interactions_bl_id.loc[cath_pdb_residue_interactions_bl_id._merge != "both"]) == 0)
        cath_pdb_residue_interactions_bl_id.drop(columns = "_merge", inplace = True)
        cath_pdb_residue_interactions_bl_id.to_csv(f"{args.outdir}/cath_pdb_residue_interactions_bl_id.csv.gz", index = False, compression = "gzip")

        scop_pdb_residue_interactions_bl_id = scop_pdb_residue_interactions_bl.merge(bound_ligands_to_score[["bound_ligand_name", "ligand_entity_id"]].rename(columns = {"ligand_entity_id": "ligand_uniqueID"}), on = "bound_ligand_name", how = "left", indicator = True)
        assert(len(scop_pdb_residue_interactions_bl_id.loc[scop_pdb_residue_interactions_bl_id._merge != "both"]) == 0)
        scop_pdb_residue_interactions_bl_id.drop(columns = "_merge", inplace = True)
        scop_pdb_residue_interactions_bl_id.to_csv(f"{args.outdir}/scop_pdb_residue_interactions_bl_id.csv.gz", index = False, compression = "gzip")

        pfam_pdb_residue_interactions_bl_id = pfam_pdb_residue_interactions_bl.merge(bound_ligands_to_score[["bound_ligand_name", "ligand_entity_id"]].rename(columns = {"ligand_entity_id": "ligand_uniqueID"}), on = "bound_ligand_name", how = "left", indicator = True)
        assert(len(pfam_pdb_residue_interactions_bl_id.loc[pfam_pdb_residue_interactions_bl_id._merge != "both"]) == 0)
        pfam_pdb_residue_interactions_bl_id.drop(columns = "_merge", inplace = True)
        pfam_pdb_residue_interactions_bl_id.to_csv(f"{args.outdir}/pfam_pdb_residue_interactions_bl_id.csv.gz", index = False, compression = "gzip")

        interpro_pdb_residue_interactions_bl_id = interpro_pdb_residue_interactions_bl.merge(bound_ligands_to_score[["bound_ligand_name", "ligand_entity_id"]].rename(columns = {"ligand_entity_id": "ligand_uniqueID"}), on = "bound_ligand_name", how = "left", indicator = True)
        assert(len(interpro_pdb_residue_interactions_bl_id.loc[interpro_pdb_residue_interactions_bl_id._merge != "both"]) == 0)
        interpro_pdb_residue_interactions_bl_id.drop(columns = "_merge", inplace = True)
        interpro_pdb_residue_interactions_bl_id.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_bl_id.csv.gz", index = False, compression = "gzip")

        bound_ligands_to_score.rename(columns = {"bound_ligand_name": "bl_name"}, inplace = True)
        bound_ligands_to_score.to_pickle(f"{args.outdir}/bound_ligands_to_score.pkl")
    else:
        bound_ligands_to_score = pd.read_pickle(f"{args.outdir}/bound_ligands_to_score.pkl")
        cath_pdb_residue_interactions_bl_id = pd.read_csv(f"{args.outdir}/cath_pdb_residue_interactions_bl_id.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        scop_pdb_residue_interactions_bl_id = pd.read_csv(f"{args.outdir}/scop_pdb_residue_interactions_bl_id.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        pfam_pdb_residue_interactions_bl_id = pd.read_csv(f"{args.outdir}/pfam_pdb_residue_interactions_bl_id.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        interpro_pdb_residue_interactions_bl_id = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_bl_id.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        console.print(f"Loaded bound ligands to score from file {args.outdir}/bound_ligands_to_score.pkl and bound ligand residue interaction ID files")

    if not os.path.exists(f"{args.outdir}/bound_entities_to_score.pkl"):
        bound_entities_to_score = pd.concat([bound_sugars_to_score[["ligand_entity_id", "bl_name", "description", "descriptor", "ec_list"]], bound_ligands_to_score])
        assert(bound_entities_to_score.ligand_entity_id.value_counts().max() == 1)
        bound_entities_to_score.to_pickle(f"{args.outdir}/bound_entities_to_score.pkl")

        cath_pdb_residue_interactions = pd.concat([cath_pdb_residue_interactions_bl_id, cath_pdb_residue_interactions_bs_index])
        cath_pdb_residue_interactions.to_csv(f"{args.outdir}/cath_pdb_residue_interactions.csv.gz", index = False, compression = "gzip")
        scop_pdb_residue_interactions = pd.concat([scop_pdb_residue_interactions_bl_id, scop_pdb_residue_interactions_bs_index])
        scop_pdb_residue_interactions.to_csv(f"{args.outdir}/scop_pdb_residue_interactions.csv.gz", index = False, compression = "gzip")
        pfam_pdb_residue_interactions = pd.concat([pfam_pdb_residue_interactions_bl_id, pfam_pdb_residue_interactions_bs_index])
        pfam_pdb_residue_interactions.to_csv(f"{args.outdir}/pfam_pdb_residue_interactions.csv.gz", index = False, compression = "gzip")
        interpro_pdb_residue_interactions = pd.concat([interpro_pdb_residue_interactions_bl_id, interpro_pdb_residue_interactions_bs_index])
        interpro_pdb_residue_interactions.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions.csv.gz", index = False, compression = "gzip")


    else:
        bound_entities_to_score = pd.read_pickle(f"{args.outdir}/bound_entities_to_score.pkl")
        cath_pdb_residue_interactions = pd.read_csv(f"{args.outdir}/cath_pdb_residue_interactions.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        scop_pdb_residue_interactions = pd.read_csv(f"{args.outdir}/scop_pdb_residue_interactions.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        pfam_pdb_residue_interactions = pd.read_csv(f"{args.outdir}/pfam_pdb_residue_interactions.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        interpro_pdb_residue_interactions = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        console.print(f"Loaded bound entities to score from file {args.outdir}/bound_entities_to_score.pkl and residue interaction files")

if __name__ == "__main__":
    main()