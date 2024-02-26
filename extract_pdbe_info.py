#!/usr/bin/env python
from pdbecif.mmcif_io import CifFileReader
from Bio.ExPASy import Enzyme as EEnzyme
from pathlib import Path
import requests
import os
import numpy as np
from utils import get_terminal_record, get_csdb_from_glycoct, get_glycoct_from_wurcs, get_smiles_from_csdb
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
from neo4j import __version__ as neo4j_version,  GraphDatabase

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
    ec_list = ec_records_df.ID.unique()
    
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

def complete_unmatched_domains(df, class_codes, fold_codes, superfamily_codes):
    df = df.merge(class_codes, left_on = "scop_class_id", right_on = "cl_id", how = "left", indicator = True)
    df["cl_description_x"] = df["cl_description_x"].fillna(df["cl_description_y"])
    df["cl_id_x"] = df["cl_id_x"].fillna(df["scop_class_id"])
    df.rename(columns = {"cl_id_x" : "cl_id", "cl_description_x": "cl_description"}, inplace = True)
    df.drop(columns = ["_merge", "cl_description_y", "cl_id_y"], inplace = True)
    df = df.merge(fold_codes, left_on = "scop_fold_id", right_on = "cf_id", how = "left", indicator = True)
    df["cf_description_x"] = df["cf_description_x"].fillna(df["cf_description_y"])
    df["cf_id_x"] = df["cf_id_x"].fillna(df["scop_fold_id"])
    df.rename(columns = {"cf_id_x" : "cf_id", "cf_description_x": "cf_description"}, inplace = True)
    df.drop(columns = [ "_merge", "cf_description_y", "cf_id_y"], inplace = True)
    df = df.merge(superfamily_codes, left_on = "scop_superfamily_id", right_on = "sf_id", how = "left", indicator = True)
    df["sf_description_x"] = df["sf_description_x"].fillna(df["sf_description_y"])
    df["sf_id_x"] = df["sf_id_x"].fillna(df["scop_superfamily_id"])
    df.rename(columns = {"sf_id_x" : "sf_id", "sf_description_x": "sf_description"}, inplace = True)
    df.drop(columns = ["_merge", "sf_description_y", "sf_id_y"], inplace = True)
    return df

#class is adapted from https://towardsdatascience.com/neo4j-cypher-python-7a919a372be7
class Neo4jConnection:
    
    def __init__(self, uri, user, pwd):
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            print("Failed to create the driver:", e)
        
    def close(self):
        if self.__driver is not None:
            self.__driver.close()
        
    def query(self, query, db=None, **kwargs):
        assert self.__driver is not None, "Driver not initialized!"
        session = None
        response = None
        try: 
            session = self.__driver.session(database=db) if db is not None else self.__driver.session() 
            response = list(session.run(query, **kwargs))
        except Exception as e:
            print("Query failed:", e)
        finally: 
            if session is not None:
                session.close()
        return response

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
        "CATH" : pdbe_graph_queries["cath_bl_query"],
        "SCOP" : pdbe_graph_queries["scop_bl_query"],
        "InterProDomain" : pdbe_graph_queries["interpro_d_bl_query"],
        "InterProFamily" : pdbe_graph_queries["interpro_f_bl_query"],
        "InterProHomologousSuperfamily" : pdbe_graph_queries["interpro_h_bl_query"]}

    bs_queries = {
        "CATH" : pdbe_graph_queries["cath_sugar_query"],
        "SCOP" : pdbe_graph_queries["scop_sugar_query"],
        "InterProDomain" : pdbe_graph_queries["interpro_d_sugar_query"],
        "InterProFamily" : pdbe_graph_queries["interpro_f_sugar_query"],
        "InterProHomologousSuperfamily" : pdbe_graph_queries["interpro_h_sugar_query"]}

    console.print("Connecting to neo4j")
    conn = Neo4jConnection(uri="bolt://localhost:7687", user="neo4j", pwd="yTJutYQ$$d%!9h")
    console.print("Connected to neo4j")
    console.print("Generating EC record dataframe")
    
    with Progress() as progress:
            
        with open(f"{args.enzyme_dat_file}") as handle:
            ec_records = EEnzyme.parse(handle)
            ec_records_list = []
            for record in ec_records: 
                ec_record_series = pd.Series(record)
                ec_records_list.append(ec_record_series)


        ec_records_df = pd.DataFrame(ec_records_list)
        ec_records_df["TRANSFER"] = ec_records_df.apply(lambda x: get_terminal_record(x["ID"], x, ec_records_df), axis = 1)
        ec_records_df["TRANSFER"] = ec_records_df["TRANSFER"].fillna(ec_records_df.ID)

        sifts_chains = pd.read_csv(f"{args.sifts_ec_mapping}", sep = "\t", comment="#")
        sifts_chains_ec = sifts_chains.loc[sifts_chains.EC_NUMBER != "?"]

        sifts_chains_ec = get_updated_enzyme_records(sifts_chains_ec, ec_records_df, ec_col = "EC_NUMBER")
        sifts_chains_ec.rename(columns = {"EC_NUMBER": "protein_entity_ec", "ACCESSION" : "uniprot_accession"}, inplace = True)
        

        all_pdbs = sifts_chains_ec.PDB.unique()
        total_rows = len(all_pdbs)

        bl_results = {}
        bs_results = {}

        scop_domains_info = pd.read_csv(f"{args.scop_domains_info_file}", sep = "\t", comment = "#", header = None, names = ["scop_id", "pdb_id", "description", "sccs", "domain_sunid", "ancestor_sunid"])
        scop_id_levels = ["cl_id", "cf_id", "sf_id", "fa_id", "dm_id", "sp_id", "px_id"]
        scop_domains_info[scop_id_levels] = scop_domains_info.ancestor_sunid.str.split(",", expand = True)
        scop_descriptions = pd.read_csv(f"{args.scop_descriptions_file}", sep = "\t", comment = "#" , header = None, names = ["level_sunid", "level", "level_sccs", "level_sid", "level_description"])

        for column in scop_id_levels:
            scop_domains_info = clean_and_merge_scop_col(scop_domains_info, column, scop_descriptions)
        
        scop_domains_info.drop(columns = ["pdb_id"], inplace = True)

        class_codes = scop_domains_info[["cl_id", "cl_description"]].drop_duplicates()
        fold_codes = scop_domains_info[["cf_id", "cf_description"]].drop_duplicates()
        superfamily_codes = scop_domains_info[["sf_id", "sf_description"]].drop_duplicates()

        if not os.path.exists(f"{args.outdir}/bl_results.pkl"):

            for db, query in bl_queries.items():
                task = progress.add_task(f"[cyan]Processing {db} bound ligands...", total=total_rows)
                results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
                    futures = [executor.submit(process_row, conn = conn, progress = progress, task = task, pdb = pdb, query = query) for pdb in all_pdbs]
                    for future in concurrent.futures.as_completed(futures):
                        if future.result() != None:
                            results.extend(future.result())

                # Convert results to DataFrame
                result_df = pd.DataFrame([dict(_) for _ in results])
                result_df = result_df.merge(sifts_chains_ec, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left", indicator = True)
                assert(len(result_df.loc[result_df._merge != "both"]) == 0)
                result_df.drop(columns = "_merge", inplace = True)
                if db == "SCOP":
                    result_df = result_df.merge(scop_domains_info, how = "left", on = "scop_id", indicator = True)
                    scop_bl_domains_matched = result_df.loc[result_df._merge == "both"].copy().drop(columns = ["_merge"])
                    scop_bl_domains_unmatched = result_df.loc[result_df._merge != "both"].copy().drop(columns = ["_merge"])
                    scop_bl_domains_unmatched = complete_unmatched_domains(scop_bl_domains_unmatched, class_codes, fold_codes, superfamily_codes)
                    result_df = pd.concat([scop_bl_domains_matched, scop_bl_domains_unmatched])
                bl_results[db] = result_df
            del(results)
            with open(f"{args.outdir}/bl_results.pkl", 'wb') as f:
                pickle.dump(bl_results, f)

        else:
            #result_df = pd.read_csv(f"{args.outdir}/{db}_pdb_residue_interactions_distinct_bl.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
            with open(f"{args.outdir}/bl_results.pkl", 'rb') as f:
                bl_results = pickle.load(f)
            console.print(f"Loaded bound ligand results from file {args.outdir}/bl_results.pkl")

        if not os.path.exists(f"{args.outdir}/bs_results.pkl"):
            for db, query in bs_queries.items():
                task = progress.add_task(f"[cyan]Processing {db} bound sugars...", total=total_rows)
                results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
                    futures = [executor.submit(process_row, conn = conn, progress = progress, task = task, pdb = pdb, query = query) for pdb in all_pdbs]
                    for future in concurrent.futures.as_completed(futures):
                        if future.result() != None:
                            results.extend(future.result())

                # Convert results to DataFrame
                result_df = pd.DataFrame([dict(_) for _ in results])
                result_df = result_df.merge(sifts_chains_ec, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left", indicator = True)
                assert(len(result_df.loc[result_df._merge != "both"]) == 0)
                result_df.drop(columns = "_merge", inplace = True)
                result_df["ligand_entity_id_numerical"] = result_df["ligand_entity_id_numerical"].astype(int)
                if db == "SCOP":
                    result_df = result_df.merge(scop_domains_info, how = "left", on = "scop_id", indicator = True)
                    scop_bl_domains_matched = result_df.loc[result_df._merge == "both"].copy().drop(columns = ["_merge"])
                    scop_bl_domains_unmatched = result_df.loc[result_df._merge != "both"].copy().drop(columns = ["_merge"])
                    scop_bl_domains_unmatched = complete_unmatched_domains(scop_bl_domains_unmatched, class_codes, fold_codes, superfamily_codes)
                    result_df = pd.concat([scop_bl_domains_matched, scop_bl_domains_unmatched])
                bs_results[db] = result_df
            del(results)
            with open(f"{args.outdir}/bs_results.pkl", 'wb') as f:
                pickle.dump(bs_results, f)
        else:
            with open(f"{args.outdir}/bs_results.pkl", 'rb') as f:
                bs_results = pickle.load(f)
            console.print(f"Loaded bound sugar results from file {args.outdir}/bs_results.pkl")

    cath_pdb_residue_interactions_bl = bl_results["CATH"]
    scop_pdb_residue_interactions_bl = bl_results["SCOP"]
    interpro_pdb_residue_interactions_bl = pd.concat([bl_results["InterProDomain"], bl_results["InterProFamily"], bl_results["InterProHomologousSuperfamily"]])
    
    #check what we are using this for again
    bound_molecules_ligands = pd.concat([cath_pdb_residue_interactions_bl[["bm_uniqids", "bound_ligand_id"]].drop_duplicates(), scop_pdb_residue_interactions_bl[["bm_uniqids", "bound_ligand_id"]].drop_duplicates(),
            interpro_pdb_residue_interactions_bl[["bm_uniqids", "bound_ligand_id"]].drop_duplicates()])
    
    bound_molecules_ligands.to_csv(f"{args.outdir}/bound_molecules_ligands.csv.gz", index = False, compression = "gzip")
    
    cath_pdb_residue_interactions_bs = bs_results["CATH"]
    scop_pdb_residue_interactions_bs = bs_results["SCOP"]
    interpro_pdb_residue_interactions_bs = pd.concat([bs_results["InterProDomain"], bs_results["InterProFamily"], bs_results["InterProHomologousSuperfamily"]])
     
    if not os.path.exists(f"{args.outdir}/bound_molecules_sugars_wurcs.csv.gz"):
        bound_molecules_sugars = pd.concat([
                cath_pdb_residue_interactions_bs[["pdb_id", "bm_uniqids", "ligand_entity_id", "uniqueID", "description", "ligand_entity_id_numerical", "protein_entity_ec", "ec_list"]].drop_duplicates() , 
                scop_pdb_residue_interactions_bs[["pdb_id", "bm_uniqids", "ligand_entity_id", "uniqueID", "description", "ligand_entity_id_numerical", "protein_entity_ec", "ec_list"]].drop_duplicates() ,
                interpro_pdb_residue_interactions_bs[["pdb_id", "bm_uniqids", "ligand_entity_id", "uniqueID", "description", "ligand_entity_id_numerical", "protein_entity_ec", "ec_list"]].drop_duplicates()
                ])

        cif_ids = bound_molecules_sugars.pdb_id.unique()
        Path(f"{args.outdir}/sugar_cifs").mkdir(parents=True, exist_ok=True)
        task = progress.add_task(f"[cyan]Downloading bound sugars CIF files...", total=len(cif_ids))
        for cif_id in cif_ids:
            cif_path = f"{args.outdir}/sugar_cifs/{cif_id}_updated.cif"
            if not os.path.exists(cif_path):
                response = requests.get(f'http://www.ebi.ac.uk/pdbe/entry-files/download/{cif_id}_updated.cif')

                with open(cif_path, 'wb') as fp:
                    fp.write(response.content)
            progress.update(task, advance=1)


        branched_entity_list = []
        reader = CifFileReader()
        task = progress.add_task(f"[cyan]Extracting WURCS from bound sugars CIF files...", total=len(cif_ids))
        for cif_id in cif_ids:
            cif_path = f'{args.outdir}/sugar_cifs/{cif_id}_updated.cif'
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
    
    if not os.path.exists(f"{args.outdir}/bound_sugars_to_score.pkl"):
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

        interpro_pre_shape = len(interpro_pdb_residue_interactions_bs)
        interpro_pdb_residue_interactions_bs_index = interpro_pdb_residue_interactions_bs.merge(bound_molecules_sugars_smiles[["ligand_index", "uniqueID", "descriptor"]].drop_duplicates().rename(columns = {"ligand_index": "ligand_uniqueID"}), on = "uniqueID", how = "left", indicator = True, validate = "many_to_one")
        assert(len(interpro_pdb_residue_interactions_bs_index) == interpro_pre_shape)
        assert(len(interpro_pdb_residue_interactions_bs_index.loc[interpro_pdb_residue_interactions_bs_index._merge != "both"]) == 0)
        interpro_pdb_residue_interactions_bs_index.drop(columns = ["_merge"], inplace = True)
        interpro_pdb_residue_interactions_bs_index.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_bs_index.csv.gz", index = False, compression = "gzip")
        
    else:
        bound_molecules_sugars_smiles = pd.read_pickle(f"{args.outdir}/bound_molecules_sugars_smiles.pkl")
        bound_sugars_to_score = pd.read_pickle(f"{args.outdir}/bound_sugars_to_score.pkl")
        console.print(f"Loaded bound molecules sugars smiles from file {args.outdir}/bound_molecules_sugars_smiles.csv.gz and bound sugars to score from {args.outdir}/bound_sugars_to_score.pkl")
    
    
    if not os.path.exists(f"{args.outdir}/bound_ligands_to_score.pkl"):
        bound_ligands_to_score = pd.concat([cath_pdb_residue_interactions_bl[["description", "bound_ligand_name", "descriptor", "ec_list"]],
                                            scop_pdb_residue_interactions_bl[["description", "bound_ligand_name", "descriptor", "ec_list"]],
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

        interpro_pdb_residue_interactions_bl_id = interpro_pdb_residue_interactions_bl.merge(bound_ligands_to_score[["bound_ligand_name", "ligand_entity_id"]].rename(columns = {"ligand_entity_id": "ligand_uniqueID"}), on = "bound_ligand_name", how = "left", indicator = True)
        assert(len(interpro_pdb_residue_interactions_bl_id.loc[interpro_pdb_residue_interactions_bl_id._merge != "both"]) == 0)
        interpro_pdb_residue_interactions_bl_id.drop(columns = "_merge", inplace = True)
        interpro_pdb_residue_interactions_bl_id.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_bl_id.csv.gz", index = False, compression = "gzip")

        bound_ligands_to_score.rename(columns = {"bound_ligand_name": "bl_name"}, inplace = True)
        bound_ligands_to_score.to_pickle(f"{args.outdir}/bound_ligands_to_score.pkl")
    else:
        bound_ligands_to_score = pd.read_pickle(f"{args.outdir}/bound_ligands_to_score.pkl")
        console.print(f"Loaded bound ligands to score from file {args.outdir}/bound_ligands_to_score.pkl")

    if not os.path.exists(f"{args.outdir}/bound_entities_to_score.pkl"):
        bound_entities_to_score = pd.concat([bound_sugars_to_score[["ligand_entity_id", "bl_name", "description", "descriptor", "ec_list"]], bound_ligands_to_score])
        assert(bound_entities_to_score.ligand_entity_id.value_counts().max() == 1)
        bound_entities_to_score.to_pickle(f"{args.outdir}/bound_entities_to_score.pkl")
    else:
        bound_entities_to_score = pd.read_pickle(f"{args.outdir}/bound_entities_to_score.pkl")
        console.print(f"Loaded bound entities to score from file {args.outdir}/bound_entities_to_score.pkl")

if __name__ == "__main__":
    main()