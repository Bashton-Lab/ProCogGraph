#!/usr/bin/env python
import pandas as pd
import json
import requests
import os
import ast
import argparse
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
from Bio.ExPASy import Enzyme as EEnzyme
import argparse
from pathlib import Path
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from utils import get_terminal_record, get_csdb_from_glycoct
from bs4 import BeautifulSoup
from urllib.parse import quote

def get_kegg_enzymes(ec_list, enzyme_string_file = None):
    def extract_reaction(enzyme_record):
        reaction_list = enzyme_record.reaction
        rn_numbers = []
        for reaction_str in reaction_list:
                reaction_ids = re.findall(r'\[RN:(.*)\]', reaction_str)
                for reaction in reaction_ids:
                    reactions = set(reaction.split())
                    rn_numbers.extend(reactions)
        if len(rn_numbers) > 0:
            return rn_numbers
        else:
            return np.nan
        
    def extract_compound_codes(text):
        pattern = r'\[CPD:([^\]]+)\]'
        matches = re.findall(pattern, text)
        if len(matches) > 0: 
            return matches
        else:
            return np.nan    
    if not enzyme_string_file:
        response = requests.get(f'https://rest.kegg.jp/get/{"+".join(ec_list)}')
        if response.status_code == 200:
            response_string = response.text
        else:
            response_string = ""
            for ec in ec_list:
                enzyme_dict[ec] = {"entry" : ec, "error" : f"KEGG API returned status code {response.status_code}"}
            return enzyme_dict, response_string
    else:
        with open(enzyme_string_file, "r") as file:
            response_string = file.read()
    enzyme_dict = {}
    
    enzyme_list = list(Enzyme.parse(io.StringIO(response_string)))
    enzyme_dict = {item.entry : {"entry": item.entry, 
                                    "error" : np.nan, 
                                    "matched_name" : item.name[0], 
                                    "reaction_text" : item.reaction, 
                                    "EC_substrate_codes" : extract_compound_codes(','.join(item.substrate)) if isinstance(item.substrate, list) else [],
                                    "EC_product_codes" : extract_compound_codes(','.join(item.product)) if isinstance(item.product, list) else [],
                                    "EC_dbxrefs": item.dblinks, 
                                    "EC_reactions" : extract_reaction(item)} for item in enzyme_list}
                
    if set(enzyme_dict.keys()) != set(ec_list):
        missing_ecs = list(set(ec_list) - set(enzyme_dict.keys()))
        for missing_ec in missing_ecs:
            enzyme_dict[missing_ec] = {"entry" : missing_ec, "error" : "KEGG API did not return result"}
    
    return enzyme_dict, response_string

def extract_secondary_id(identifier, database_list, current_db = ""):
    #takes the identifier, and the database list
    #if the identifier is a glycan , preference is to search for the structure in the glycosmos dataset.
    #else, preference is [CHEBI,
    
    pattern = r'^G'
    match = re.match(pattern, identifier)
    
    database = None
    identifier = None
    for item in database_list:
        if database not in ["GlyTouCan", "ChEBI"]:
            if match and item[0] == "GlyTouCan" and item[0] != current_db:
                database = item[0]
                identifier = item[1][0]
            else:
                if not match and item[0] == "ChEBI" and item[0] != current_db:
                    database = item[0]
                    identifier = item[1][0]
                elif not match and item[0] == "PubChem":
                    database = item[0]
                    identifier = item[1][0]
    
    return database, identifier

def get_kegg_reactions(chunk, reactions_string_file = None):
    kegg_reaction_dictionary = {}
    if not reactions_string_file:
        response = requests.get(f'https://rest.kegg.jp/get/{ "+".join(chunk)}')
        if response.status_code == 200:
            responses_string = response.text
        else:
            for reaction in chunk:
                kegg_reaction_dictionary[reaction] = {"reaction_id": reaction, "reaction_definition" : f"KEGG Reaction not found (Error: {response.status_code}"}
            responses_string = ""
            return kegg_reaction_dictionary, responses_string
    else:
        with open(reactions_string_file, "r") as file:
            responses_string = file.read()
        
    responses_string_split = responses_string.split("///")
    for item in responses_string_split:
        record = item.split("\n")
        for line in record:
            if line.startswith('ENTRY'):
                entry = line.split()[1]
            elif line.startswith('DEFINITION'):
                definition = line.split(maxsplit = 1)[1:]
            elif line.startswith('EQUATION'):
                equation = line.split(maxsplit = 1)[1:][0]
            elif line.startswith('ENZYME'):
                enzyme = line.split(maxsplit = 1)[1:]
        substrates, products = equation.split(' <=> ')
        substrate_codes = re.findall(r'\b[C|G]\d+\b', substrates)
        product_codes = re.findall(r'\b[C|G]\d+\b', products)
        kegg_reaction_dictionary[entry] = {"reaction_id": entry, 
                                            "reaction_definition" : definition, 
                                            "reaction_equation": equation, 
                                            "reaction_substrates": substrates, 
                                            "reaction_products" : products, 
                                            "reaction_substrate_codes" : substrate_codes, 
                                            "reaction_product_codes" : product_codes}
        for reaction in chunk:
            if reaction not in kegg_reaction_dictionary.keys():
                kegg_reaction_dictionary[reaction] = {"reaction_id": reaction, "reaction_definition" : f"KEGG Reaction not found (Error: {response.status_code}"}
    
    return kegg_reaction_dictionary, responses_string


#pubchem listkeys may speed this up? https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial#section=Dealing-with-Lists-of-Identifiers
#may subsequently want to change property to InChI or IsomericSMILES down the line. 
def pubchem_cid_to_descriptor(compound_list, chunk_size = 200):
    df_list = []
    for i in range(0, len(compound_list), chunk_size):
        chunk = compound_list[i:i+chunk_size]
        compound_string = ",".join(chunk)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{compound_string}/property/CanonicalSMILES/JSON"
        # Fetch JSON data from the URL
        response = requests.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            # Load JSON data
            data = response.json().get("PropertyTable", {}).get("Properties", [])
            df = pd.DataFrame(data)
            df_list.append(df)
        else:
            # If the request was unsuccessful, print an error message
            print("Failed to fetch data. Status code:", response.status_code)
        
    # Convert to DataFrame
    df = pd.concat(df_list)

    return df

def get_kegg_compound_record(kegg_id, glycan = False): 
    compound_record = requests.get(f'https://rest.kegg.jp/get/{kegg_id}')
    if compound_record.status_code == 200:
        compound_record_object = list(Enzyme.parse(io.StringIO(compound_record.text)))[0]
        dblinks = compound_record_object.dblinks
        if not glycan:
            name = compound_record_object.name[0]
        else:
            name = None
    else:
        compound_record_object = None
        dblinks = None
        name = None
        
    compound_dict = {"compound_id" : kegg_id, "compound_name": name, "dbxrefs": dblinks, "KEGG_compound_record" : compound_record_object}
    return compound_dict

def get_kegg_compound_smiles(kegg_id):
    response = requests.get(f'https://rest.kegg.jp/get/{kegg_id}/mol')
    if response.status_code == 200:
        compound_split = response.text.split("> <ENTRY>\n")
        molblock = compound_split[0]
        if molblock == "":   
            canonical_smiles = np.nan
        else:
            canonical_smiles = canon_smiles(Chem.MolFromMolBlock(molblock))
    else:
        canonical_smiles = np.nan
    return canonical_smiles

def get_smiles_from_csdb(csdb_linear):
    if csdb_linear is np.nan:
        return np.nan
    else:
        response = requests.get(f"http://csdb.glycoscience.ru/database/core/convert_api.php?csdb={quote(csdb_linear)}&format=smiles")
        mol = np.nan
        smiles = np.nan
        if response.status_code == 200:
            html = response.text
            soup = BeautifulSoup(html, 'html.parser')
            for a in soup.find_all("a"):
                title = a.get('title')
                if title == "Find this structure in ChemSpider":
                    smiles = a.contents[0].strip()
                    break
        else:
            smiles = np.nan   
        return smiles
    
def get_gtc_info(gtcids):
    glytoucan_df_list = []
    for gtcid in gtcids:
        url = 'https://api.glycosmos.org/sparqlist/gtcid2seqs'
        data = {'gtcid': gtcid}
        nested_dict = {}

        response = requests.post(url, data=data)

        if response.status_code == 200:
            json_result = json.loads(response.text)  # Display the response content
            # Iterate through each dictionary in the list
            for item in json_result:
                id_val = item['id']
                wurcs_val = item.get('wurcs')
                glycoct_val = item.get('glycoct')

                # Check if the ID exists in the nested dictionary
                if id_val not in nested_dict:
                    # Create a new entry for the ID
                    nested_dict[id_val] = {'wurcs': wurcs_val, 'glycoct': glycoct_val}
                else:
                    # Update the existing ID entry with non-None values
                    if wurcs_val is not None:
                        nested_dict[id_val]['wurcs'] = wurcs_val
                    if glycoct_val is not None:
                        nested_dict[id_val]['glycoct'] = glycoct_val
            df_dict = pd.DataFrame(nested_dict).T
            glytoucan_df_list.append(df_dict)
        else:
            print(f"Request failed with status code: {response.status_code}")
            
    gtc_df = pd.concat(glytoucan_df_list)

    gtc_df.reset_index(inplace = True)
    return gtc_df

def canon_smiles(x):
        try:
            return Chem.MolToSmiles(x)
        except:
            return np.nan
        
from rdkit.Chem import PandasTools

# Function to calculate descriptors
def calculate_descriptors(id, mol, calculator):
    return (id,) + calculator.CalcDescriptors(mol)

def main():    
    """
    This script is designed to take the enzyme.dat file from EXPASY, and extract the EC numbers, 
    and then use the KEGG API to extract the associated EC records, reaction records and compound records.
    The script then searches databases: KEGG, RHEA, ChEBI, PubChem and GlyTouCan for the compound codes, and
    extracts the SMILES strings for each compound. Where possible, duplicate structures are aggregated
    based on their canonical SMILES into a single biological ligand record (with a list of ligand_db identifiers).
    The script then uses the RDKit to calculate molecular descriptors for each ligand, and saves the results
    in a pkl format.
    """

    parser = argparse.ArgumentParser(description='Process EC information.')
    parser.add_argument('--ec_dat', type=str, help='Path to enzyme.dat file from EXPASY')
    parser.add_argument('--pubchem', type=str, help='Path to pubchem_substance_id_mapping file')
    parser.add_argument('--chebi', type=str, help='Path to chebi_kegg_file')
    parser.add_argument('--rhea_mapping', type=str, help='Path to rhea-directions.tsv file')
    parser.add_argument('--rhea_reactions', type=str, help='Path to rhea-reaction-smiles.tsv file')
    parser.add_argument('--rhea2ec', type=str, help='Path to rhea2ec.tsv file')
    parser.add_argument('--outdir', type=str, help='Path to output directory')

    args = parser.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    
    if os.path.exists(f"{args.outdir}/biological_ligands.pkl"):
        print("Biological ligands file already exists. Exiting.")
        exit(0)

    with open(args.ec_dat) as handle:
        ec_records = EEnzyme.parse(handle)
        ec_records_list = []
        for record in ec_records: 
            ec_record_series = pd.Series(record)
            ec_records_list.append(ec_record_series)
            
    ec_records_df = pd.DataFrame(ec_records_list)
    ec_records_df["TRANSFER"] = ec_records_df.apply(lambda x: get_terminal_record(x["ID"], x, ec_records_df), axis = 1)
    ec_records_df["TRANSFER"] = ec_records_df["TRANSFER"].fillna(ec_records_df.ID)

    ec_list = ec_records_df.TRANSFER.unique()

    if not os.path.exists(f"{args.outdir}/kegg_enzyme_df.pkl"):
        print("Getting KEGG Enzyme records")
        n=10 #chunk size
        enzyme_records = {}
        for i in range(0, len(ec_list), n):
            print(f"Processing chunk {i} of {len(ec_list)}")
            chunk = ec_list[i:i + n]
            enzyme_dict = get_kegg_enzymes(chunk)
            enzyme_records.update(enzyme_dict)
            time.sleep(1)

        kegg_enzyme_df = pd.DataFrame(enzyme_records).T
        kegg_enzyme_df = kegg_enzyme_df.explode("EC_reactions")
        kegg_enzyme_df.drop(columns = ['reaction_text', 'EC_dbxrefs'], inplace = True)
        kegg_enzyme_df.to_pickle(f"{args.outdir}/kegg_enzyme_df.pkl")
        print("KEGG Enzyme records saved")
    else:
        kegg_enzyme_df = pd.read_pickle(f"{args.outdir}/kegg_enzyme_df.pkl")
        print("KEGG Enzyme records loaded from file")
    reactions = kegg_enzyme_df.EC_reactions.dropna().unique()

    if not os.path.exists(f"{args.outdir}/kegg_reaction_df.pkl"):
        print("Getting KEGG Reaction records")
        kegg_reaction_dictionary = {}
        n=10 #chunk size
        for i in range(0, len(reactions), n):
            chunk = reactions[i:i + n]
            reaction_dictionary = get_kegg_reactions(chunk)
            kegg_reaction_dictionary.update(reaction_dictionary)
            time.sleep(1)

        kegg_reaction_df = pd.DataFrame(kegg_reaction_dictionary).T
        assert(len(kegg_reaction_df.loc[kegg_reaction_df.reaction_definition.str.len() > 1]) == 0) #what is this assertion for? 
        kegg_reaction_df = kegg_reaction_df.explode("reaction_definition")
        kegg_reaction_df["reaction_definition"] = kegg_reaction_df.reaction_id + ": " + kegg_reaction_df.reaction_definition
        kegg_reaction_df["reaction_equation"] = kegg_reaction_df.reaction_id + ": " + kegg_reaction_df.reaction_equation
        kegg_reaction_df.to_pickle(f"{args.outdir}/kegg_reaction_df.pkl")
        print("KEGG Reaction records saved")
    else:
        kegg_reaction_df = pd.read_pickle(f"{args.outdir}/kegg_reaction_df.pkl")    
        print("KEGG Reaction records loaded from file")
    ##need to have some sort of verification here that the reaction df is based on the same reactions that the enzyme df is based on.

    if not os.path.exists(f"{args.outdir}/kegg_reaction_enzyme_df.pkl"):
        print("Merging KEGG Enzyme and Reaction records")
        kegg_reaction_enzyme_df = kegg_enzyme_df.merge(kegg_reaction_df, left_on = "EC_reactions", right_on = "reaction_id", how = "left")
        kegg_reaction_enzyme_df["reaction_substrate_codes"] = kegg_reaction_enzyme_df["reaction_substrate_codes"].apply(lambda d: d if isinstance(d, list) else [])
        kegg_reaction_enzyme_df["reaction_product_codes"] = kegg_reaction_enzyme_df["reaction_product_codes"].apply(lambda d: d if isinstance(d, list) else [])
        kegg_reaction_enzyme_df["EC_substrate_codes"] = kegg_reaction_enzyme_df["EC_substrate_codes"].apply(lambda d: d if isinstance(d, list) else [])
        kegg_reaction_enzyme_df["EC_product_codes"] = kegg_reaction_enzyme_df["EC_product_codes"].apply(lambda d: d if isinstance(d, list) else [])

        kegg_reaction_enzyme_df = kegg_reaction_enzyme_df.fillna("").groupby("entry").agg({"entry" : set, "error" : set, "matched_name" : set, "EC_substrate_codes": sum,"EC_product_codes": sum, "reaction_substrate_codes" : sum, "reaction_product_codes": sum, "EC_reactions":set, "reaction_id" : set, "reaction_definition": set, "reaction_equation" : set})
        kegg_reaction_enzyme_df["entities"] = kegg_reaction_enzyme_df["EC_substrate_codes"] + kegg_reaction_enzyme_df["EC_product_codes"] + kegg_reaction_enzyme_df["reaction_substrate_codes"] + kegg_reaction_enzyme_df["reaction_product_codes"]
        kegg_reaction_enzyme_df["entities"] = kegg_reaction_enzyme_df["entities"].apply(lambda x: ",".join(list(set(x))))
        #convert substrate/product codes to a comma separated string instead of list (for easier merging later)
        kegg_reaction_enzyme_df["EC_substrate_codes"] = kegg_reaction_enzyme_df["EC_substrate_codes"].apply(lambda x: ",".join(x))
        kegg_reaction_enzyme_df["EC_product_codes"] = kegg_reaction_enzyme_df["EC_product_codes"].apply(lambda x: ",".join(x))
        kegg_reaction_enzyme_df["reaction_substrate_codes"] = kegg_reaction_enzyme_df["reaction_substrate_codes"].apply(lambda x: ",".join(x))
        kegg_reaction_enzyme_df["reaction_product_codes"] = kegg_reaction_enzyme_df["reaction_product_codes"].apply(lambda x: ",".join(x))
        kegg_reaction_enzyme_df["reaction_definition"] = kegg_reaction_enzyme_df["reaction_definition"].apply(lambda x: ",".join(x))
        kegg_reaction_enzyme_df["reaction_equation"] = kegg_reaction_enzyme_df["reaction_equation"].apply(lambda x: ",".join(x))

        def unpack_sets(row):
            for col in row.index:
                if isinstance(row[col], set) and len(row[col]) == 1:
                    row[col] = next(iter(row[col]))
                elif isinstance(row[col], set) and len(row[col]) > 1:
                    row[col] = list(row[col])
                elif isinstance(row[col], set) and len(row[col]) == 0:
                    row[col] = np.nan
            return row

        kegg_reaction_enzyme_df = kegg_reaction_enzyme_df.apply(unpack_sets, axis=1)
        kegg_reaction_enzyme_df.to_pickle(f"{args.outdir}/kegg_reaction_enzyme_df.pkl")
        print("KEGG Enzyme and Reaction records merged and saved")
    else:
        kegg_reaction_enzyme_df = pd.read_pickle(f"{args.outdir}/kegg_reaction_enzyme_df.pkl")
        print("Merged KEGG Enzyme and Reaction records loaded from file")

    kegg_reaction_enzyme_df["entities"] = kegg_reaction_enzyme_df.entities.str.split(",", expand = False)

    #use the coompound codes to filter the pubchem data.
    kegg_reaction_enzyme_df_exploded = kegg_reaction_enzyme_df.explode("entities")
    compound_codes = kegg_reaction_enzyme_df_exploded.entities.dropna().unique()

    #here we should try and get the kegg compounds from the compound codes, where they are non glycans (deal with glycans later)
    if not os.path.exists(f"{args.outdir}/kegg_compounds_df.pkl"):
        print("Getting KEGG Compound records")
        compound_codes = kegg_reaction_enzyme_df_exploded.entities.dropna().unique()
        kegg_compounds_df = pd.DataFrame([get_kegg_compound_record(code) for code in compound_codes if code.startswith("C")])
        kegg_compounds_df["canonical_smiles"] = kegg_compounds_df["compound_id"].apply(lambda x: get_kegg_compound_smiles(x))
        kegg_compounds_df.to_pickle(f"{args.outdir}/kegg_compounds_df.pkl")
        print("KEGG Compound records saved")
    else:
        kegg_compounds_df = pd.read_pickle(f"{args.outdir}/kegg_compounds_df.pkl")
        print("KEGG Compound records loaded from file")

    if not os.path.exists(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_kegg.pkl"):
        print("Merging KEGG Compound records with Enzyme and Reaction records")
        kegg_reaction_enzyme_df_exploded_kegg = kegg_reaction_enzyme_df_exploded.merge(kegg_compounds_df, left_on = "entities", right_on = "compound_id", how = "inner")
        kegg_reaction_enzyme_df_exploded_kegg = kegg_reaction_enzyme_df_exploded_kegg.dropna(subset = ["canonical_smiles"])
        kegg_reaction_enzyme_df_exploded_kegg["ligand_db"] = "KEGG:" + kegg_reaction_enzyme_df_exploded_kegg["entities"].astype("str")
        PandasTools.AddMoleculeColumnToFrame(kegg_reaction_enzyme_df_exploded_kegg, smilesCol='canonical_smiles')
        kegg_reaction_enzyme_df_exploded_kegg.to_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_kegg.pkl")
        print("KEGG Compound records merged with Enzyme and Reaction records and saved")
    else:
        kegg_reaction_enzyme_df_exploded_kegg = pd.read_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_kegg.pkl")
        print("Merged KEGG Compound records with Enzyme and Reaction records loaded from file")


    #the rhea reactions can have identical molecules on left and right hand side, so should drop duplicates for final df.
    if not os.path.exists(f"{args.outdir}/kegg_enzyme_df_rhea.pkl"):
        print("Getting RHEA records")
        rhea_mapping = pd.read_csv(f"{args.rhea_mapping}", sep = "\t")
        rhea_reactions = pd.read_csv(f"{args.rhea_reactions}", sep = "\t", header = None, names = ["REACTION_ID", "REACTION"])
        rhea2ec = pd.read_csv(f"{args.rhea2ec}", sep = "\t")

        rhea_mapping_merged = rhea_mapping.merge(rhea2ec, left_on = "RHEA_ID_MASTER", right_on = "MASTER_ID", how = "right")

        #we use the left-to-right mapping of the reaction for identifying smiles strings for the compounds
        rhea_mapping_merged_reactions = rhea_mapping_merged.merge(rhea_reactions, left_on = "RHEA_ID_LR", right_on = "REACTION_ID", how = "left")
        rhea_mapping_merged_reactions = rhea_mapping_merged_reactions.loc[rhea_mapping_merged_reactions.REACTION.isna() == False].copy()
        kegg_enzyme_df_rhea = kegg_enzyme_df.merge(rhea_mapping_merged_reactions, left_on = "entry", right_on = "ID", how = "inner")
        
        kegg_enzyme_df_rhea = kegg_enzyme_df_rhea.loc[kegg_enzyme_df_rhea.REACTION.isna() == False].copy()
        kegg_enzyme_df_rhea["REACTANT"] = kegg_enzyme_df_rhea["REACTION"].str.split(r">>|\.")
        kegg_enzyme_df_rhea = kegg_enzyme_df_rhea.explode("REACTANT")
        kegg_enzyme_df_rhea = kegg_enzyme_df_rhea.drop_duplicates(subset = ["entry", "REACTANT"])
        kegg_enzyme_df_rhea.reset_index(drop = True, inplace = True)
        PandasTools.AddMoleculeColumnToFrame(kegg_enzyme_df_rhea, smilesCol='REACTANT')

        kegg_enzyme_df_rhea = kegg_enzyme_df_rhea.loc[kegg_enzyme_df_rhea.ROMol.isna() == False].reset_index(drop = True)
        kegg_enzyme_df_rhea["ligand_db"] = "RHEA:" + kegg_enzyme_df_rhea["RHEA_ID"].astype("str")
        kegg_enzyme_df_rhea.to_pickle(f"{args.outdir}/kegg_enzyme_df_rhea.pkl")
        print("RHEA records saved")
    else:
        kegg_enzyme_df_rhea = pd.read_pickle(f"{args.outdir}/kegg_enzyme_df_rhea.pkl")
        print("RHEA records loaded from file")
    ### ChEBI KEGG Mapping
    if not os.path.exists(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_chebi.pkl"):
        print("Getting ChEBI records")
        chebi_kegg_compounds = pd.read_csv(args.chebi, sep = "\t", index_col = False)
        chebi_kegg_compounds = chebi_kegg_compounds.rename(columns = {"ID" : "ChEBI_ID", "NAME" : "ChEBI_NAME"})
        chebi_kegg_compounds_smiles = chebi_kegg_compounds.loc[chebi_kegg_compounds.SMILES.isna() == False]

        ### Merging ChEBI with enzyme df

        kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded.merge(chebi_kegg_compounds_smiles, left_on = "entities", right_on = "KEGG COMPOUND ACCESSION", how = "inner", indicator = True)
        kegg_reaction_enzyme_df_exploded_chebi.drop(columns = "_merge", inplace = True)
        kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded_chebi.loc[kegg_reaction_enzyme_df_exploded_chebi.SMILES.isna() == False].copy()

        PandasTools.AddMoleculeColumnToFrame(kegg_reaction_enzyme_df_exploded_chebi, smilesCol='SMILES')
        kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded_chebi.loc[kegg_reaction_enzyme_df_exploded_chebi.ROMol.isna() == False].reset_index(drop = True)

        kegg_reaction_enzyme_df_exploded_chebi["ligand_db"] = kegg_reaction_enzyme_df_exploded_chebi["ChEBI_ID"].astype("str")
        kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded_chebi.loc[kegg_reaction_enzyme_df_exploded_chebi.ROMol.isna() == False].reset_index(drop = True)
        kegg_reaction_enzyme_df_exploded_chebi.to_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_chebi.pkl")
        print("ChEBI records saved")
    else:
        kegg_reaction_enzyme_df_exploded_chebi = pd.read_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_chebi.pkl")    
        print("ChEBI records loaded from file")

    if not os.path.exists(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_pubchem.pkl"):
        print("Getting PubChem records")
        ### PubChem KEGG Mapping from here : https://pubchem.ncbi.nlm.nih.gov/source/
        with open(args.pubchem, 'r') as file:
            lines = file.readlines()
            
        data = []
        for line in lines:
            record = {}
            line = line.replace(": ", ":")
            parts = line.strip().split()  # Split line by space
            for part in parts:
                key, value = part.split(':')  # Split by colon and space to get key-value pairs
                record[key] = value
            data.append(record)

        kegg_pubchem_mapping = pd.DataFrame(data)

        cids = kegg_pubchem_mapping.loc[(kegg_pubchem_mapping.CID.isna() == False) &
                                        (kegg_pubchem_mapping.KEGG.isin(compound_codes))].CID.unique()

        pubchem_smiles = pubchem_cid_to_descriptor(cids)

        kegg_pubchem_mapping["CID"] = kegg_pubchem_mapping["CID"].fillna(0).astype(int)
        pubchem_smiles["CID"] = pubchem_smiles["CID"].astype(int)

        kegg_pubchem_mapping_smiles = kegg_pubchem_mapping.merge(pubchem_smiles, on = "CID", how = "left", indicator = True)
        kegg_pubchem_mapping_smiles_filtered = kegg_pubchem_mapping_smiles.loc[kegg_pubchem_mapping_smiles.CanonicalSMILES.isna() == False].reset_index(drop = True).copy()
        kegg_pubchem_mapping_smiles_filtered.drop(columns = "_merge", inplace = True)

        ### Merging PubChem with enzyme df
        kegg_reaction_enzyme_df_exploded_pubchem = kegg_reaction_enzyme_df_exploded.merge(kegg_pubchem_mapping_smiles_filtered, left_on = "entities", right_on = "KEGG", how = "inner")
        PandasTools.AddMoleculeColumnToFrame(kegg_reaction_enzyme_df_exploded_pubchem, smilesCol='CanonicalSMILES')

        kegg_reaction_enzyme_df_exploded_pubchem = kegg_reaction_enzyme_df_exploded_pubchem.loc[kegg_reaction_enzyme_df_exploded_pubchem.ROMol.isna() == False].reset_index(drop = True)
        kegg_reaction_enzyme_df_exploded_pubchem["ligand_db"] = "Pubchem:" + kegg_reaction_enzyme_df_exploded_pubchem["CID"].astype("str")
        kegg_reaction_enzyme_df_exploded_pubchem.to_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_pubchem.pkl")
    else:
        kegg_reaction_enzyme_df_exploded_pubchem = pd.read_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_pubchem.pkl")


    if not os.path.exists(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_gtc.pkl"):
        ### Getting list of glytoucan glycans to query
        glycans = kegg_reaction_enzyme_df_exploded.loc[(kegg_reaction_enzyme_df_exploded.entities.isna() == False) 
                                                        & (kegg_reaction_enzyme_df_exploded.entities.str.startswith("G")), "entities"].values.tolist()

        ### Getting KEGG compound records for glycans (to get xref to glytoucan)
        glycan_compounds = []

        for glycan in glycans:
            compound = get_kegg_compound_record(glycan, glycan = True)
            glycan_compounds.append(compound)
            
        glycan_compounds_df = pd.DataFrame(glycan_compounds)
        glycan_compounds_df[["source", "secondary_id"]] = glycan_compounds_df.apply(lambda x: extract_secondary_id(x["compound_id"], x["dbxrefs"]), axis = 1, result_type = "expand").values
        glytoucan_ids = glycan_compounds_df.loc[(glycan_compounds_df.secondary_id.isna() == False) &
                                        (glycan_compounds_df.source == "GlyTouCan"), "secondary_id"].unique()

        glytoucan_df = get_gtc_info(glytoucan_ids)

        glycan_compounds_df_merged = glycan_compounds_df.merge(glytoucan_df, left_on = "secondary_id", right_on = "index", how = "left")

        #now need to convert this through into an smiles representation.
        #convert to csdb linear and then smiles
        glycan_compounds_df_merged["csdb_linear"] = glycan_compounds_df_merged.glycoct.apply(lambda x: get_csdb_from_glycoct(x))
        glycan_compounds_df_merged["smiles"] = glycan_compounds_df_merged.apply(lambda x: get_smiles_from_csdb(x["csdb_linear"]), axis = 1)
        glycan_compounds_df_merged = glycan_compounds_df_merged.loc[glycan_compounds_df_merged.smiles.isna() == False]
        glycan_compounds_df_merged["name"] = glycan_compounds_df_merged["secondary_id"]
        kegg_reaction_enzyme_df_exploded_gtc = kegg_reaction_enzyme_df_exploded.merge(glycan_compounds_df_merged, left_on="entities", right_on="compound_id", how = "inner")
        PandasTools.AddMoleculeColumnToFrame(kegg_reaction_enzyme_df_exploded_gtc, smilesCol='smiles')

        kegg_reaction_enzyme_df_exploded_gtc = kegg_reaction_enzyme_df_exploded_gtc.loc[kegg_reaction_enzyme_df_exploded_gtc.ROMol.isna() == False].reset_index(drop = True)
        kegg_reaction_enzyme_df_exploded_gtc["ligand_db"] = "GlyTouCan:" + kegg_reaction_enzyme_df_exploded_gtc["secondary_id"].astype("str")
        kegg_reaction_enzyme_df_exploded_gtc.to_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_gtc.pkl")
        print("GlyTouCan records saved")
    else:
        print("GlyTouCan records loaded from file")
        kegg_reaction_enzyme_df_exploded_gtc = pd.read_pickle(f"{args.outdir}/kegg_reaction_enzyme_df_exploded_gtc.pkl")
        
    #append together the dataframes with various representations of biological ligands into a long form dataframe. Should have the EC number, EC name
    #the kegg compound entry ID, the db identifier for the biological molecule, the unique identifier for the biological ligand molecule, and the descriptor (mol/smiles/inchi)
    #this will be used to generate the parity scores for the biological ligands. 

    #How do we make a unique identifier for the molecules?
    #try to also get a name for the compounds.

    #maybe we want to get the name from the kegg entities records for all instances, and then merge that to get a name
    #need to clean up the dataframes here to remove : "_merge" before import to neo4j file script

    if not os.path.exists(f"{args.outdir}/biological_ligands_df.pkl"):
        descriptors = ['BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 
                'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'ExactMolWt', 'FpDensityMorgan1', 'FpDensityMorgan2', 
                'FpDensityMorgan3', 'FractionCSP3', 'HallKierAlpha', 'HeavyAtomCount', 'HeavyAtomMolWt', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'MaxAbsEStateIndex', 'MaxAbsPartialCharge', 
                'MaxEStateIndex', 'MaxPartialCharge', 'MinAbsEStateIndex', 'MinAbsPartialCharge', 'MinEStateIndex', 'MinPartialCharge', 'MolLogP', 'MolMR', 'MolWt', 'NHOHCount', 'NOCount', 
                'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 
                'NumHDonors', 'NumHeteroatoms', 'NumRadicalElectrons', 'NumRotatableBonds', 'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 'NumSaturatedRings', 'NumValenceElectrons', 
                'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 
                'PEOE_VSA9', 'RingCount', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 
                'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'TPSA', 'VSA_EState1', 'VSA_EState10', 
                'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 
                'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 
                'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 
                'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine', 'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether', 'fr_furan', 
                'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 
                'fr_methoxy', 'fr_morpholine', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 
                'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 'fr_sulfonamd', 
                'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea', 'qed']
        #kegg_enzyme_df_rhea[["entry", "ROMol", "ligand_db"]],  removign rhea from concat - if adding back 
        # create molecular descriptor calculator
        mol_descriptor_calculator = MolecularDescriptorCalculator(descriptors)
        biological_ligands_df = pd.concat([kegg_reaction_enzyme_df_exploded_kegg[["entry", "compound_name", "compound_id", "ROMol", "ligand_db"]], 
                                           kegg_reaction_enzyme_df_exploded_chebi[["entry", "ChEBI_NAME", "KEGG COMPOUND ACCESSION", "ROMol", "ligand_db"]].rename(columns = {"KEGG COMPOUND ACCESSION" : "compound_id"}), 
                                           kegg_reaction_enzyme_df_exploded_pubchem[["entry", "KEGG", "ROMol", "ligand_db"]].rename(columns = {"KEGG" : "compound_id"}),
                                           kegg_reaction_enzyme_df_exploded_gtc[["entry", "name", "compound_id","ROMol", "ligand_db"]].rename(columns = {"name": "compound_name"})])
        biological_ligands_df = biological_ligands_df.reset_index()
        

        biological_ligands_df["compound_name"] = biological_ligands_df["compound_name"].fillna(biological_ligands_df["ChEBI_NAME"]).fillna(biological_ligands_df["compound_id"])
        biological_ligands_df.drop(columns = ["ChEBI_NAME"], inplace = True)
        biological_ligands_df = biological_ligands_df.groupby(["entry", "compound_id"], dropna = False).agg({"ROMol": "first", "compound_name" : "first", "ligand_db": set}).reset_index()
        biological_ligands_df["ligand_db"] = biological_ligands_df["ligand_db"].str.join("|")
        biological_ligands_df["canonical_smiles"] = biological_ligands_df["ROMol"].map(lambda x: canon_smiles(x) if isinstance(x,Chem.rdchem.Mol) else np.nan)
        biological_ligands_df = biological_ligands_df.groupby(["entry", "canonical_smiles", "compound_id", "compound_name"], dropna = False).agg({"ligand_db": set, "ROMol" : "first"}).reset_index()
        biological_ligands_df = biological_ligands_df.drop_duplicates(subset = ["entry", "canonical_smiles"])
        biological_ligands_df["uniqueID"] = biological_ligands_df.groupby('canonical_smiles').ngroup()
        biological_ligands_df["ligand_db"] = biological_ligands_df["ligand_db"].str.join("|")
        #biological_ligands_mol_descriptors = biological_ligands_df['ROMol'].apply(calculate_descriptors)
        #biological_ligands_mol_descriptors_df = pd.DataFrame(list(biological_ligands_mol_descriptors), columns=descriptors)
        #biological_ligands_df = pd.concat([biological_ligands_df, biological_ligands_mol_descriptors_df], axis=1)
        biological_ligands_df.to_pickle(f"{args.outdir}/biological_ligands_df.pkl")
    else:
        biological_ligands_df = pd.read_pickle(f"{args.outdir}/biological_ligands_df.pkl")

if __name__ == "__main__":
    main()
    
#glygen??