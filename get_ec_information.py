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

from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

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

def canon_smiles(x):
        try:
            return Chem.MolToSmiles(x)
        except:
            return np.nan
        
from rdkit.Chem import PandasTools

# Function to calculate descriptors
def calculate_descriptors(id, mol):
    return (id,) + mol_descriptor_calculator.CalcDescriptors(mol)

parser = argparse.ArgumentParser(description='Process EC information.')
parser.add_argument('--ec_dat', type=str, help='Path to enzyme.dat file from EXPASY')
parser.add_argument('--pubchem', type=str, help='Path to pubchem_substance_id_mapping file')
parser.add_argument('--chebi', type=str, help='Path to chebi_kegg_file')
parser.add_argument('--reaction_enzyme', type=str, help='')

args = parser.parse_args()


descriptors = ['BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'ExactMolWt', 'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'FractionCSP3', 'HallKierAlpha', 'HeavyAtomCount', 'HeavyAtomMolWt', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'MaxAbsEStateIndex', 'MaxAbsPartialCharge', 'MaxEStateIndex', 'MaxPartialCharge', 'MinAbsEStateIndex', 'MinAbsPartialCharge', 'MinEStateIndex', 'MinPartialCharge', 'MolLogP', 'MolMR', 'MolWt', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRadicalElectrons', 'NumRotatableBonds', 'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 'NumSaturatedRings', 'NumValenceElectrons', 'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'RingCount', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'TPSA', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine', 'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea', 'qed']

# create molecular descriptor calculator
mol_descriptor_calculator = MolecularDescriptorCalculator(descriptors)

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
    kegg_enzyme_df = kegg_enzyme_df.explode("EC_reactions")
    kegg_enzyme_df.drop(columns = ['reaction_text', 'EC_dbxrefs'], inplace = True)
    kegg_enzyme_df.to_pickle("kegg_enzyme_df.pkl")

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
    assert(len(kegg_reaction_df.loc[kegg_reaction_df.reaction_definition.str.len() > 1]) == 0)
    kegg_reaction_df = kegg_reaction_df.explode("reaction_definition")
    kegg_reaction_df["reaction_definition"] = kegg_reaction_df.reaction_id + ": " + kegg_reaction_df.reaction_definition
    kegg_reaction_df["reaction_equation"] = kegg_reaction_df.reaction_id + ": " + kegg_reaction_df.reaction_equation
    kegg_reaction_df.to_pickle("kegg_reaction_df.pkl")

try:
    kegg_reaction_enzyme_df = pd.read_pickle(args.reaction_enzyme)
except:
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
    kegg_reaction_enzyme_df.to_pickle("kegg_reaction_enzyme_df.pkl")

kegg_reaction_enzyme_df["entities"] = kegg_reaction_enzyme_df.entities.str.split(",", expand = False)

compound_codes = kegg_reaction_enzyme_df.entities.explode().dropna().unique()
kegg_reaction_enzyme_df_exploded = kegg_reaction_enzyme_df.explode("entities")
compound_codes = kegg_reaction_enzyme_df_exploded.entities.dropna().unique()

chebi_compounds = pd.read_csv("kegg_information/chebi_compounds.tsv", sep = "\t")
#filter to remove anything where name is not present
chebi_compounds_filtered = chebi_compounds.loc[chebi_compounds.NAME.isna() == False]

# def get_compound_id_by_name(data_frame, name_to_match, compounds = "", relations = "", structures = "",search = False):
#     print(f"working on: {name_to_match}")
#     name_to_match_singular = re.sub(r's$', '', name_to_match)
#     name_to_match_a = f"a {name_to_match}"
#     name_to_match_an = f"an {name_to_match}"
    
#     match_list = [name_to_match, name_to_match_singular,name_to_match_a, name_to_match_an]
#     matching_row = data_frame.loc[data_frame['NAME'].str.lower().isin([x.lower() for x in match_list])]
    
#     if not matching_row.empty:
#         compound_id = matching_row['ID'].astype("int").iloc[0]
#         if search:
#             compound_structure = find_chebi_structure(compound_id, compounds, relations, structures)
#             return compound_structure
#         else:
#             return compound_id
#     else:
#         if search:
#             return [["", "no structure identified", ""]]
#         else:
#             return None
# print(chebi_compounds)

# print(kegg_reaction_enzyme_df_exploded.loc[(kegg_reaction_enzyme_df_exploded.entities.isna()) & (kegg_reaction_enzyme_df_exploded.EC_reactions.isna()) & 
#                                 (kegg_reaction_enzyme_df_exploded.EC_substrate_codes == "") & 
#                                 (kegg_reaction_enzyme_df_exploded.EC_product_codes == "") & 
#                                 (kegg_reaction_enzyme_df_exploded.matched_name != "Deleted entry"),  ["entry","reaction_text"]].drop_duplicates(["entry"]).dropna(subset = ["reaction_text"]).copy())

#should add the ROmol info to the molecule side before merging. Reduces the amount of duplications potentially.
#check the merging and how it is working. 
#we probably want an inner join? keep only the intersection of keys.

#the rhea reactions can have identical molecules on left and right hand side, so should drop duplicates for final df.
try:
    kegg_enzyme_df_rhea = pd.read_pickle("kegg_enzyme_df_rhea.pkl")
except:
    rhea_mapping = pd.read_csv("rhea-directions.tsv", sep = "\t")
    rhea_reactions = pd.read_csv("rhea-reaction-smiles.tsv", sep = "\t", header = None, names = ["REACTION_ID", "REACTION"])
    rhea2ec = pd.read_csv("rhea2ec.tsv", sep = "\t")

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

    kegg_enzyme_df_rhea["canonical_smiles"] = kegg_enzyme_df_rhea["ROMol"].map(lambda x: canon_smiles(x) if isinstance(x,Chem.rdchem.Mol) else np.nan)

    # Calculate descriptors for each molecule in the 'ROMol' column
    #rhea_mol_descriptors = kegg_enzyme_df_rhea['ROMol'].apply(calculate_descriptors)
    #rhea_mol_descriptors_df = pd.DataFrame(list(rhea_mol_descriptors), columns=descriptors)

    # Merge the descriptors DataFrame with the 'test' DataFrame
    #kegg_enzyme_df_rhea = pd.concat([kegg_enzyme_df_rhea, rhea_mol_descriptors_df], axis=1)
    kegg_enzyme_df_rhea["ligand_db"] = "RHEA"
    kegg_enzyme_df_rhea = kegg_enzyme_df_rhea.loc[kegg_enzyme_df_rhea.ROMol.isna() == False].reset_index(drop = True)
    kegg_enzyme_df_rhea.to_pickle("kegg_enzyme_df_rhea.pkl")

### ChEBI KEGG Mapping
try:
    kegg_reaction_enzyme_df_exploded_chebi = pd.read_pickle("kegg_reaction_enzyme_df_exploded_chebi.pkl")
except:
    chebi_kegg_compounds = pd.read_csv(args.chebi, sep = "\t", index_col = False)
    chebi_kegg_compounds = chebi_kegg_compounds.rename(columns = {"ID" : "ChEBI_ID", "NAME" : "ChEBI_NAME"})
    chebi_kegg_compounds_smiles = chebi_kegg_compounds.loc[chebi_kegg_compounds.SMILES.isna() == False]

    ### Merging ChEBI with enzyme df

    kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded.merge(chebi_kegg_compounds_smiles, left_on = "entities", right_on = "KEGG COMPOUND ACCESSION", how = "inner", indicator = True)
    kegg_reaction_enzyme_df_exploded_chebi.drop(columns = "_merge", inplace = True)
    kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded_chebi.loc[kegg_reaction_enzyme_df_exploded_chebi.SMILES.isna() == False].copy()

    PandasTools.AddMoleculeColumnToFrame(kegg_reaction_enzyme_df_exploded_chebi, smilesCol='SMILES')
    kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded_chebi.loc[kegg_reaction_enzyme_df_exploded_chebi.ROMol.isna() == False].reset_index(drop = True)
    kegg_reaction_enzyme_df_exploded_chebi["canonical_smiles"] = kegg_reaction_enzyme_df_exploded_chebi["ROMol"].map(lambda x: canon_smiles(x) if isinstance(x,Chem.rdchem.Mol) else np.nan)

    # Calculate descriptors for each molecule in the 'ROMol' column
    #chebi_mol_descriptors = kegg_reaction_enzyme_df_exploded_chebi['ROMol'].apply(calculate_descriptors)
    #chebi_mol_descriptors_df = pd.DataFrame(list(chebi_mol_descriptors), columns=descriptors)
    # Merge the descriptors DataFrame with the 'test' DataFrame
    #kegg_reaction_enzyme_df_exploded_chebi = pd.concat([kegg_reaction_enzyme_df_exploded_chebi, chebi_mol_descriptors_df], axis=1)
    kegg_reaction_enzyme_df_exploded_chebi["ligand_db"] = "ChEBI"
    kegg_reaction_enzyme_df_exploded_chebi = kegg_reaction_enzyme_df_exploded_chebi.loc[kegg_reaction_enzyme_df_exploded_chebi.ROMol.isna() == False].reset_index(drop = True)
    kegg_reaction_enzyme_df_exploded_chebi.to_pickle("kegg_reaction_enzyme_df_exploded_chebi.pkl")
    

try:
    kegg_reaction_enzyme_df_exploded_pubchem = pd.read_pickle("kegg_reaction_enzyme_df_exploded_pubchem.pkl")
except:
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
    kegg_reaction_enzyme_df_exploded_pubchem.rename(columns={"CanonicalSMILES" : "canonical_smiles"}, inplace = True)
    PandasTools.AddMoleculeColumnToFrame(kegg_reaction_enzyme_df_exploded_pubchem, smilesCol='canonical_smiles')

    # Calculate descriptors for each molecule in the 'ROMol' column
    #pubchem_mol_descriptors = kegg_reaction_enzyme_df_exploded_pubchem['ROMol'].apply(calculate_descriptors)
    #pubchem_mol_descriptors_df = pd.DataFrame(list(pubchem_mol_descriptors), columns=descriptors)

    # Merge the descriptors DataFrame with the 'test' DataFrame
    #kegg_reaction_enzyme_df_exploded_pubchem = pd.concat([kegg_reaction_enzyme_df_exploded_pubchem, pubchem_mol_descriptors_df], axis=1)
    kegg_reaction_enzyme_df_exploded_pubchem["ligand_db"] = "PubChem"
    kegg_reaction_enzyme_df_exploded_pubchem.rename(columns = {"CanonicalSMILES" : "canonical_smiles"}, inplace = True)
    kegg_reaction_enzyme_df_exploded_pubchem = kegg_reaction_enzyme_df_exploded_pubchem.loc[kegg_reaction_enzyme_df_exploded_pubchem.ROMol.isna() == False].reset_index(drop = True)
    kegg_reaction_enzyme_df_exploded_pubchem.to_pickle("kegg_reaction_enzyme_df_exploded_pubchem.pkl")

try:
    kegg_reaction_enzyme_df_exploded_gtc = pd.read_pickle("kegg_reaction_enzyme_df_exploded_gtc.pkl")
    
except:
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
    glycan_compounds_df_merged = glycan_compounds_df_merged.loc[glycan_compounds_df_merged.smiles.isna() == False]
    kegg_reaction_enzyme_df_exploded_gtc = kegg_reaction_enzyme_df_exploded.merge(glycan_compounds_df_merged, left_on="entities", right_on="compound_id", how = "inner")
    PandasTools.AddMoleculeColumnToFrame(kegg_reaction_enzyme_df_exploded_gtc, smilesCol='smiles')

    kegg_reaction_enzyme_df_exploded_gtc["canonical_smiles"] = kegg_reaction_enzyme_df_exploded_gtc["ROMol"].map(lambda x: canon_smiles(x) if isinstance(x,Chem.rdchem.Mol) else np.nan)

    # Calculate descriptors for each molecule in the 'ROMol' column
    #gtc_mol_descriptors = kegg_reaction_enzyme_df_exploded_gtc['ROMol'].apply(calculate_descriptors)
    #gtc_mol_descriptors_df = pd.DataFrame(list(gtc_mol_descriptors), columns=descriptors)

    # Merge the descriptors DataFrame with the 'test' DataFrame
    #kegg_reaction_enzyme_df_exploded_gtc = pd.concat([kegg_reaction_enzyme_df_exploded_gtc, gtc_mol_descriptors_df], axis=1)
    kegg_reaction_enzyme_df_exploded_gtc = kegg_reaction_enzyme_df_exploded_gtc.drop_duplicates(subset = ["entry", "canonical_smiles"])
    kegg_reaction_enzyme_df_exploded_gtc["ligand_db"] = "GlyTouCan"
    kegg_reaction_enzyme_df_exploded_gtc = kegg_reaction_enzyme_df_exploded_gtc.loc[kegg_reaction_enzyme_df_exploded_gtc.ROMol.isna() == False].reset_index(drop = True)
    kegg_reaction_enzyme_df_exploded_gtc.to_pickle("kegg_reaction_enzyme_df_exploded_gtc.pkl")

#append together the dataframes with various representations of biological ligands into a long form dataframe. Should have the EC number, EC name
#the kegg compound entry ID, the db identifier for the biological molecule, the unique identifier for the biological ligand molecule, and the descriptor (mol/smiles/inchi)
#this will be used to generate the parity scores for the biological ligands. 

#How do we make a unique identifier for the molecules?

#try to also get a name for the compounds.



#biological_ligand_df
biological_ligands_df = pd.concat([kegg_reaction_enzyme_df_exploded_chebi[["entry", "canonical_smiles", "ROMol", "ligand_db"]], kegg_reaction_enzyme_df_exploded_pubchem[["entry", "canonical_smiles", "ROMol", "ligand_db"]], kegg_enzyme_df_rhea[["entry", "canonical_smiles", "ROMol", "ligand_db"]], kegg_reaction_enzyme_df_exploded_gtc[["entry", "canonical_smiles", "ROMol", "ligand_db"]]])
biological_ligands_df = biological_ligands_df.reset_index()

biological_ligands_df_databases = biological_ligands_df[["entry", "canonical_smiles", "ligand_db"]].groupby(["entry", "canonical_smiles"]).agg(list).reset_index()

biological_ligands_df = biological_ligands_df.drop_duplicates(subset = ["entry", "canonical_smiles"]).drop(columns = "ligand_db")
biological_ligands_df = biological_ligands_df.merge(biological_ligands_df_databases, on = ["entry", "canonical_smiles"], how = "left", indicator = True)
print(biological_ligands_df.loc[biological_ligands_df._merge != "both"])
assert(len(biological_ligands_df.loc[biological_ligands_df._merge != "both"]) == 0)
biological_ligands_df["uniqueID"] = biological_ligands_df.groupby('canonical_smiles').ngroup()
biological_ligands_df["ligand_db"] = biological_ligands_df["ligand_db"].str.join("|")
print(biological_ligands_df)

biological_ligands_descriptors = biological_ligands_df.apply(lambda x: calculate_descriptors(x["uniqueID"], x["ROMol"]), axis = 1)
print(biological_ligands_descriptors)
biological_ligands_descriptors_df = pd.DataFrame(list(biological_ligands_descriptors), columns=["uniqueID"] + descriptors)
biological_ligands_df = pd.concat([biological_ligands_df, biological_ligands_descriptors_df], axis=1)

biological_ligands_df.to_pickle("biological_ligands_df.pkl")
#should be 70297