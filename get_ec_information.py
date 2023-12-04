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

def get_kegg_enzymes(ec_list):
    def extract_reaction(enzyme_record):
        reaction_list = enzyme_record.reaction
        rn_numbers = []
        for reaction_str in reaction_list:
            rn_numbers.extend(re.findall(r'\[RN:(R\d+)\]', reaction_str))
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
    
    response = requests.get(f'https://rest.kegg.jp/get/{"+".join(ec_list)}')
    enzyme_dict = {}
    if response.status_code == 200:
        enzyme_list = list(Enzyme.parse(io.StringIO(response.text)))
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
            
    else:
        for ec in ec_list:
            enzyme_dict[ec] = {"entry" : ec, "error" : f"KEGG API returned status code {response.status_code}"}
    
    return enzyme_dict

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

def get_kegg_reactions(chunk):
    kegg_reaction_dictionary = {}
    response = requests.get(f'https://rest.kegg.jp/get/{ "+".join(chunk)}')
    if response.status_code == 200:
        responses_string = response.text
        responses_string = responses_string.split("///")
        for item in responses_string:
            record = item.split("\n")
            for line in record:
                if line.startswith('ENTRY'):
                    entry = line.split()[1]
                elif line.startswith('DEFINITION'):
                    definition = line.split(maxsplit = 1)[1:]
                elif line.startswith('EQUATION'):
                    equation = line.split(maxsplit = 1)[1:][0]

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
    else:
        for reaction in chunk:
            kegg_reaction_dictionary[reaction] = {"reaction_id": reaction, "reaction_definition" : f"KEGG Reaction not found (Error: {response.status_code}"}
    return kegg_reaction_dictionary

def pubchem_cid_to_inchi(compound_list, chunk_size = 200):
    df_list = []
    for i in range(0, len(compound_list), chunk_size):
        chunk = compound_list[i:i+chunk_size]
        compound_string = ",".join(chunk)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{compound_string}/property/InChI/JSON"
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

def get_kegg_compound_record(kegg_id):
    def extract_dblinks(record):
        dblinks = record.dblinks
        return dblinks
    
    compound_record = requests.get(f'https://rest.kegg.jp/get/{kegg_id}')
    if compound_record.status_code == 200:
        compound_record_object = list(Enzyme.parse(io.StringIO(compound_record.text)))[0]
        dblinks = extract_dblinks(compound_record_object)
    else:
        compound_record_object = None
        dblinks = None
        
    compound_dict = {"compound_id" : kegg_id, "dbxrefs": dblinks, "KEGG_compound_record" : compound_record_object}
    return compound_dict

def get_kegg_compound_inchi(kegg_id):
    response = requests.get(f'https://rest.kegg.jp/get/{kegg_id}/mol')
    if response.status_code == 200:
        compound_split = response.text.split("> <ENTRY>\n")
        molblock = compound_split[0]
        if molblock == "":   
            molblock = None
    else:
        molblock = None

#a limitation of this function is that it selects only one representative EC number for transferred entries
def get_terminal_record(entry, row, df):
    entry = row.ID
    pattern = r"[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+"
    while row.DE.startswith("Transferred entry: "):
        transfers = re.findall(pattern, row.DE)
        #when updating, if multiple possible transfers, selecting only first.
        row = df.loc[df.ID == transfers[0]].iloc[0]
    return row.ID


def get_csdb_from_glycoct(glycoct):
    if glycoct is np.nan:
        return np.nan
    else:
        url = "http://csdb.glycoscience.ru/database/core/convert_api.php"
        data = {"glycoct":glycoct}
        headers = {'Content-Type': 'application/json'}
        response = requests.get(f"{url}?glycoct={quote(glycoct)}")
        if response.status_code == 200:
            response_data = response.text.replace("<pre>", "")  # Remove the <pre> tag
            lines = response_data.split("\n")  # Split into lines
            csdb_linear = np.nan
            for line in lines:
                if line.startswith("CSDB Linear:"):
                    csdb_linear = line.replace("CSDB Linear:", "").strip()  # Extract the CSDB Linear string
                    break
        else:
            csdb_linear = np.nan
        return csdb_linear

from bs4 import BeautifulSoup
from urllib.parse import quote
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

parser = argparse.ArgumentParser(description='Process EC information.')
parser.add_argument('--ec_dat', type=str, help='Path to enzyme.dat file from EXPASY')
parser.add_argument('--ec_class', type=str, help='URL to enzclass.txt file')
parser.add_argument('--pubchem', type=str, help='Path to pubchem_substance_id_mapping file')
parser.add_argument('--chebi', type=str, help='Path to chebi_kegg_file')
parser.add_argument('--reaction_enzyme', type=str, help='')

args = parser.parse_args()

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

try:
    kegg_enzyme_df = pd.read_pickle("kegg_enzyme_df.pkl")
except:
    n=10 #chunk size
    enzyme_records = {}
    for i in range(0, len(ec_list), n):
        print(f"Processing chunk {i} of {len(ec_list)}")
        chunk = ec_list[i:i + n]
        enzyme_dict = get_kegg_enzymes(chunk)
        enzyme_records.update(enzyme_dict)
        time.sleep(1)

    kegg_enzyme_df = pd.DataFrame(enzyme_records).T
    kegg_enzyme_df.to_pickle("kegg_enzyme_df.pkl")
kegg_enzyme_df = kegg_enzyme_df.explode("EC_reactions")

reactions = kegg_enzyme_df.EC_reactions.dropna().unique()

try:
    kegg_reaction_df = pd.read_pickle("kegg_reaction_df.pkl")
except:
    kegg_reaction_dictionary = {}
    n=10 #chunk size
    for i in range(0, len(reactions), n):
        chunk = reactions[i:i + n]
        reaction_dictionary = get_kegg_reactions(chunk)
        kegg_reaction_dictionary.update(reaction_dictionary)
        time.sleep(1)

    kegg_reaction_df = pd.DataFrame(kegg_reaction_dictionary).T
    kegg_reaction_df.to_pickle("kegg_reaction_df.pkl")


try:
    kegg_reaction_enzyme_df = pd.read_pickle(args.reaction_enzyme)
except:
    kegg_reaction_enzyme_df = kegg_enzyme_df.merge(kegg_reaction_df, left_on = "EC_reactions", right_on = "reaction_id", how = "left")
    kegg_reaction_enzyme_df["reaction_substrate_codes"] = kegg_reaction_enzyme_df["reaction_substrate_codes"].apply(lambda d: d if isinstance(d, list) else [])
    kegg_reaction_enzyme_df["reaction_product_codes"] = kegg_reaction_enzyme_df["reaction_product_codes"].apply(lambda d: d if isinstance(d, list) else [])
    kegg_reaction_enzyme_df["EC_substrate_codes"] = kegg_reaction_enzyme_df["EC_substrate_codes"].apply(lambda d: d if isinstance(d, list) else [])
    kegg_reaction_enzyme_df["EC_product_codes"] = kegg_reaction_enzyme_df["EC_product_codes"].apply(lambda d: d if isinstance(d, list) else [])

    kegg_reaction_enzyme_df["entities"] = kegg_reaction_enzyme_df["EC_substrate_codes"] + kegg_reaction_enzyme_df["EC_product_codes"] + kegg_reaction_enzyme_df["reaction_substrate_codes"] + kegg_reaction_enzyme_df["reaction_product_codes"]
    kegg_reaction_enzyme_df["entities"] = kegg_reaction_enzyme_df["entities"].apply(lambda x: list(set(x)))

    kegg_reaction_enzyme_df.to_pickle("kegg_reaction_enzyme_df.pkl")


compound_codes = kegg_reaction_enzyme_df.entities.explode().dropna().unique()
            
kegg_reaction_enzyme_df_exploded = kegg_reaction_enzyme_df.explode("entities")
compound_codes = kegg_reaction_enzyme_df_exploded.entities.dropna().unique()

### ChEBI KEGG Mapping
chebi_kegg_compounds = pd.read_csv(args.chebi, sep = "\t", index_col = False)
chebi_kegg_compounds = chebi_kegg_compounds.rename(columns = {"ID" : "ChEBI_ID", "NAME" : "ChEBI_NAME"})
chebi_kegg_compounds["CHEBI_INCHI"] = chebi_kegg_compounds["ChEBI_ID"] + "_" + chebi_kegg_compounds["INCHI"]
chebi_kegg_compounds_inchi = chebi_kegg_compounds.loc[chebi_kegg_compounds.CHEBI_INCHI.isna() == False]
chebi_kegg_compounds_inchi_grouped = chebi_kegg_compounds.groupby("KEGG COMPOUND ACCESSION").agg({"CHEBI_INCHI" : list, "GLYTOUCAN ACCESSION": set}).reset_index()
chebi_kegg_compounds_inchi_grouped["CHEBI_INCHI"] = chebi_kegg_compounds_inchi_grouped["CHEBI_INCHI"].astype(str)
### Merging ChEBI with enzyme df

kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded.merge(chebi_kegg_compounds_inchi_grouped, left_on = "entities", right_on = "KEGG COMPOUND ACCESSION", how = "left", indicator = True)
kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded_chebi.explode("GLYTOUCAN ACCESSION")
kegg_reaction_enzyme_df_exploded_chebi_exploded = kegg_reaction_enzyme_df_exploded_chebi.explode("CHEBI_INCHI")
kegg_reaction_enzyme_df_exploded_chebi_exploded.drop(columns = "_merge", inplace = True)

kegg_reaction_enzyme_df_exploded_chebi_exploded.to_pickle("kegg_reaction_enzyme_df_exploded_chebi_exploded.pkl")

### PubChem KEGG Mapping
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

pubchem_inchis = pubchem_cid_to_inchi(cids)

kegg_pubchem_mapping["CID"] = kegg_pubchem_mapping["CID"].fillna(0).astype(int)
pubchem_inchis["CID"] = pubchem_inchis["CID"].astype(int)

kegg_pubchem_mapping_inchi = kegg_pubchem_mapping.merge(pubchem_inchis, on = "CID", how = "left", indicator = True)
kegg_pubchem_mapping_inchi_filtered = kegg_pubchem_mapping_inchi.loc[kegg_pubchem_mapping_inchi.InChI.isna() == False].copy()
kegg_pubchem_mapping_inchi_filtered.drop(columns = "_merge", inplace = True)

### Merging PubChem with enzyme df
kegg_reaction_enzyme_df_exploded_pubchem = kegg_reaction_enzyme_df_exploded.merge(kegg_pubchem_mapping_inchi_filtered, left_on = "entities", right_on = "KEGG", how = "left", indicator = True)
kegg_reaction_enzyme_df_exploded_pubchem.drop(columns = "_merge", inplace = True)
kegg_reaction_enzyme_df_exploded_pubchem.to_pickle("kegg_reaction_enzyme_df_exploded_pubchem.pkl")

### Getting list of glytoucan glycans to query
glycans = kegg_reaction_enzyme_df_exploded.loc[(kegg_reaction_enzyme_df_exploded.entities.isna() == False) 
                                                & (kegg_reaction_enzyme_df_exploded.entities.str.startswith("G")), "entities"].values.tolist()

### Getting KEGG compound records for glycans (to get xref to glytoucan)
glycan_compounds = []

for glycan in glycans:
    compound = get_kegg_compound_record(glycan)
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

kegg_reaction_enzyme_df_exploded_chebi_exploded_pubchem_gtc = kegg_reaction_enzyme_df_exploded.merge(glycan_compounds_df_merged, left_on="entities", right_on="compound_id", how = "left", indicator = True)
kegg_reaction_enzyme_df_exploded_chebi_exploded_pubchem_gtc.to_pickle("biological_ligands.pkl")




