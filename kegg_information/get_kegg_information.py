import pandas as pd
import json
import requests
from jinja2 import Template
import os
from spacy import displacy
from spacy_llm.util import assemble
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

def get_kegg_compound(kegg_id):
    compound_record = requests.get(f'https://rest.kegg.jp/get/{kegg_id}')
    if compound_record.status_code == 200:
        compound_record_object = list(Enzyme.parse(io.StringIO(compound_record.text)))[0]
        dblinks = extract_dblinks(compound_record_object)
    else:
        compound_record_object = None
        dblinks = None
    if kegg_id.startswith("G"):
        response = requests.get(f'https://rest.kegg.jp/get/{kegg_id}/kcf')
        if response.status_code == 200:
            wurcs = get_wurcs_from_kcf(response.text)
            if isinstance(wurcs, str):
                glycoct = get_glycoct_from_wurcs(wurcs)
                if isinstance(glycoct, str):
                    csdb_linear = get_csdb_from_glycoct(glycoct)
                    if isinstance(csdb_linear, str):
                        mol = get_smiles_from_csdb(csdb_linear)
                    else:
                        mol = None
                else:
                    mol = None

            else:
                mol = None
        else:
            mol = None
        
    else:
        response = requests.get(f'https://rest.kegg.jp/get/{kegg_id}/mol')
        if response.status_code == 200:
            compound_split = response.text.split("> <ENTRY>\n")
            molfile = compound_split[0]
            mol = Chem.MolFromMolBlock(molfile)
            if mol == "":
                mol = None
                
        else:
            mol = None
    compound_dict = {"compound_id" : kegg_id, "mol_source": "KEGG", "mol": mol, "dbxrefs": dblinks, "KEGG_compound_record" : compound_record_object}
    return compound_dict


def extract_reaction(enzyme_record):
    reaction_list = enzyme_record.reaction
    rn_numbers = []
    for reaction_str in reaction_list:
        rn_numbers.extend(re.findall(r'\[RN:(R\d+)\]', reaction_str))
    if len(rn_numbers) > 0:
        return rn_numbers
    else:
        return np.nan

def get_kegg_enzymes(ec_list, console):
    response = requests.get(f'https://rest.kegg.jp/get/{"+".join(ec_list)}')
    if response.status_code == 200:
        enzyme_list = list(Enzyme.parse(io.StringIO(response.text)))
        enzyme_dict = {item.entry : {"pdb_entry": item.entry, "matched_entry" : item.entry, "matched_name" : item.name[0], "reaction_text" : item.reaction, "EC_substrates" : item.substrate, "EC_products" : item.product, "EC_dbxrefs": item.dblinks, "EC_reactions" : extract_reaction(item)} for item in enzyme_list}
        
        pattern = r"[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+"
        for enzyme in enzyme_dict.values():
            while ("Transferred to" in enzyme["matched_name"]) and (enzyme["matched_name"] != "Transferred to"):
                matches = re.findall(pattern, enzyme["matched_name"])
                updated_enzyme_id = next((match for match in matches if match != enzyme["pdb_entry"]), None)
                if updated_enzyme_id is not None:
                    enzyme["matched_entry"] = updated_enzyme_id
                    response = requests.get(f'https://rest.kegg.jp/get/{updated_enzyme_id}')
                    if response.status_code == 200:
                        kegg_enzyme_record = list(Enzyme.parse(io.StringIO(response.text)))[0]
                        enzyme["matched_name"] = kegg_enzyme_record.name[0]
                        enzyme["reaction_text"] = kegg_enzyme_record.reaction
                        enzyme["EC_substrates"] = kegg_enzyme_record.substrate
                        enzyme["EC_products"] = kegg_enzyme_record.product
                        enzyme["EC_dbxrefs"] = kegg_enzyme_record.dblinks
                        enzyme["EC_reactions"] = extract_reaction(kegg_enzyme_record)
                    else:
                        console.print(response.status_code)
                        enzyme["matched_entry"] = f"{updated_enzyme_id} (Error: {response.status_code})"
                else:
                    enzyme["matched_entry"] = f"Matched entry equal to original entry"
                    enzyme["matched_record"] = f"Matched entry equal to original entry"
                    
        if set(enzyme_dict.keys()) != set(ec_list):
            missing_ecs = list(set(ec_list) - set(enzyme_dict.keys()))
            for missing_ec in missing_ecs:
                enzyme_dict[missing_ec] = {"pdb_entry" : missing_ec, "matched_entry" : "KEGG API did not return result"}
            
    else:
        console.print(f"Got status code: {response.status_code}")
        for ec in ec_list:
            enzyme_dict[ec] = {"pdb_entry" : ec, "matched_entry" : f"KEGG API returned status code {response.status_code}"}
    
    return enzyme_dict

def extract_compound_codes(text):
    pattern = r'\[CPD:([^\]]+)\]'
    matches = re.findall(pattern, text)
    if len(matches) > 0: 
        return matches
    else:
        return np.nan
    

    
def extract_dblinks(record):
    dblinks = record.dblinks
    return dblinks

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

def get_kegg_reactions(chunk, console):
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
            for reaction in reactions[i:i + n]:
                if reaction not in kegg_reaction_dictionary.keys():
                    kegg_reaction_dictionary[reaction] = {"reaction_id": reaction, "reaction_definition" : f"KEGG Reaction not found (Error: {response.status_code}"}
    else:
        console.print(f"Got status code: {response.status_code}")
        for reaction in reactions[i:i + n]:
            kegg_reaction_dictionary[reaction] = {"reaction_id": reaction, "reaction_definition" : f"KEGG Reaction not found (Error: {response.status_code}"}
    return kegg_reaction_dictionary

from bs4 import BeautifulSoup

def get_smiles_from_csdb(csdb_linear):
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
        if smiles != np.nan:
            try:
                mol = Chem.MolFromSmiles(smiles)
            except:
                pass    
    return mol

from urllib.parse import quote
import json

def get_glycoct_from_wurcs(wurcs):
    url = "https://api.glycosmos.org/glycanformatconverter/2.8.2/wurcs2glycoct"
    data = {"input":wurcs}
    headers = {'Content-Type': 'application/json'}
    response = requests.post(url, headers=headers, data=json.dumps(data))
    
    if response.status_code == 200:
        response_data = response.json()
        if 'message' in response_data and response_data['message'] == 'Returned null.':
            glycoct = np.nan
        else:
            glycoct = response_data['GlycoCT']
    else:
        glycoct = np.nan
    return glycoct

def get_wurcs_from_kcf(kcf):
    url = "https://api.glycosmos.org/glycanformatconverter/2.8.2/kcf2wurcs"
    data = {"input":kcf}
    headers = {'Content-Type': 'application/json'}
    response = requests.post(url, headers=headers, data=json.dumps(data))
    
    if response.status_code == 200:
        response_data = response.json()
        if 'message' in response_data and response_data['message'] == 'Returned null.':
            wurcs = np.nan
        elif "WURCS" in response_data.keys():
            wurcs = response_data['WURCS']
        else:
            wurcs = np.nan
    else:
        wurcs = np.nan
    return wurcs

def get_csdb_from_glycoct(glycoct):
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

import pandas as pd
import xml.etree.ElementTree as ET

def parse_table_data(elem):
    data = {}
    for row in elem.findall('row'):
        row_data = {}
        for field in row.findall('field'):
            row_data[field.attrib['name']] = field.text
        data[len(data)] = row_data
    return data



def find_chebi_structure(compound_id: int, compounds, relations, structures):    
    """
    Finds the structure associated with a compound ID by calling the recursive search function.

    Args:
        compound_id (int): The compound ID to search for.
        compounds (DataFrame): DataFrame containing compound information.
        relations (DataFrame): DataFrame containing relationship information.
        structures (DataFrame): DataFrame containing structure information.

    Returns:
        dict: A dictionary of found structures mapped to their respective compound IDs,
              along with the relationship context in which they were found.

    Example:
        compound_id = 25348
        compounds_df = ...  # Your compounds DataFrame
        relations_df = ...  # Your relations DataFrame
        structures_df = ...  # Your structures DataFrame
        found_structure = find_structure(compound_id, compounds_df, relations_df, structures_df)
        print("Found structures:", found_structure)
    """
    
    def find_structure_recursive(compound_id, relationship = ""):
        """
        Recursively searches for a structure associated with a compound ID by exploring relationships.

        Args:
            compound_id (int): The compound ID to start the search from.
            relationship (str, optional): The relationship context in which the search is performed.

        Returns:
            list or None: A list of found structures mapped to their respective compound IDs,
                         along with the relationship context in which they were found. Returns None
                         if no structure is identified.

        Example:
            compound_id = 17908
            relationship = ""  # Initial relationship context
            found_structure = find_structure_recursive(compound_id, relationship)
            if found_structure is not None:
                print("Found structures:", found_structure)
            else:
                print("No structure found.")
        """
     
        # Check if the compound_id exists in the structures dataframe
        results = []

        if compound_id in structures['COMPOUND_ID'].values:
            found_structures = structures.loc[(structures['COMPOUND_ID'] == compound_id) & (structures.TYPE.isin(["mol", "InChI","SMILES"]))].copy()
            found_structures["TYPE"] = pd.Categorical(found_structures['TYPE'], ["mol", "InChI", "SMILES"])
            found_structures = found_structures.sort_values(by = "TYPE")
            if len(found_structures) > 0:
                match = found_structures.iloc[0]["STRUCTURE"]
                match_type = found_structures.iloc[0]["TYPE"]
                if match_type == "mol":
                    mol = Chem.MolFromMolBlock(match)
                elif match_type == "InChI":
                    mol = Chem.MolFromInchi(match)
                else:
                    mol = Chem.MolFromSmiles(match)
                if mol in [np.nan, None, ""]:
                    mol = None
            else:
                mol = None
            results = [[compound_id, mol, relationship]]
            return results
        # Find related compound_ids based on is_a relationships
        related_compounds = relations.loc[(relations['INIT_ID'] == compound_id) & (relations['TYPE'] == "is_a"), 'FINAL_ID'].values.tolist()
        # Recursively search for each related compound_id
        relationship = f"{relationship}:{compound_id}"
        for related_id in related_compounds:
            result = find_structure_recursive(related_id, relationship)
            if result is not None:
                results.extend(result)
        if len(results) == 0:
            return None
        else:
            return results
        
    blacklist = []#[36976]
    compound_id = int(compound_id)
    if compound_id not in blacklist:
        found_structure = find_structure_recursive(compound_id)
    else:
        found_structure = None
    if found_structure is None:
        found_structure = [[compound_id, "no structure identified", ""]]
    return(found_structure)


def get_compound_id_by_name(data_frame, name_to_match, compounds = "", relations = "", structures = "",search = False):
    print(f"working on: {name_to_match}")
    name_to_match_singular = re.sub(r's$', '', name_to_match)
    name_to_match_a = f"a {name_to_match}"
    name_to_match_an = f"an {name_to_match}"
    
    match_list = [name_to_match, name_to_match_singular,name_to_match_a, name_to_match_an]
    matching_row = data_frame.loc[data_frame['NAME'].str.lower().isin([x.lower() for x in match_list])]
    
    if not matching_row.empty:
        compound_id = matching_row['ID'].astype("int").iloc[0]
        if search:
            compound_structure = find_chebi_structure(compound_id, compounds, relations, structures)
            return compound_structure
        else:
            return compound_id
    else:
        if search:
            return [["", "no structure identified", ""]]
        else:
            return None

def get_spacy_entities(pipeline, text):
    result = nlp(text)
    entities = []

    for ent in result.ents:
        entities.append(ent.text)
    return entities

parser = argparse.ArgumentParser(description = 'TO DO')
parser.add_argument('--openai_api_key', metavar = '', type = str, 
    help = "")
parser.add_argument('--openai_api_org', metavar = '', type = str, 
    help = "")

args = parser.parse_args()


try:
    kegg_compounds_df= pd.read_pickle("glytoucan_kegg_compounds_df.pkl")
    kegg_compound_reaction_enzyme_df = pd.read_pickle("final_kegg_compound_reaction_enzyme_df.pkl")
except:
    
    # Load the XML file containing explorenz ec data
    xml_file_path = 'enzyme-data.xml'
    tree = ET.parse(xml_file_path)

    # Find all table_data elements and parse their data
    data = {}
    for table_data_elem in tree.findall(".//table_data"):
        table_name = table_data_elem.attrib['name']
        data[table_name] = parse_table_data(table_data_elem)

    # Convert the parsed data into a DataFrame
    dfs = {table_name: pd.DataFrame.from_dict(data) for table_name, data in data.items()}

    total_ec_nums = dfs["entry"].T
    total_ec_nums_list = total_ec_nums["ec_num"].dropna().unique().tolist()
    
    chains_domains = pd.read_csv("chains_domains_formatted.csv")

    ec_list = chains_domains.filled_EC_descriptor.dropna().unique().tolist()

    n=10 #chunk size

    with Progress() as progress: #transient=True

        download_enzymes = progress.add_task("[red]Downloading KEGG Enzyme Records...", total=len(range(0, len(ec_list), n)) +1)
        get_reactions = progress.add_task("[blue]Getting KEGG reaction records...", total = None)
        get_products_and_substrates = progress.add_task("[red]Downloading KEGG Compound Records...", total = None)

        while not progress.finished:
            try:
                kegg_enzyme_df = pd.read_pickle('kegg_enzyme_df.pkl')
                progress.update(download_enzymes, advance=len(range(0, len(ec_list), n)) +1)
                time.sleep(0.1)
            except:
                enzyme_records = {}
                for i in range(0, len(ec_list), n):
                    chunk = ec_list[i:i + n]
                    enzyme_dict = get_kegg_enzymes(chunk, progress.console)
                    enzyme_records.update(enzyme_dict)
                    progress.update(download_enzymes, advance=1)
                    time.sleep(1)

                kegg_enzyme_df = pd.DataFrame(enzyme_records).T
                kegg_enzyme_df['EC_substrate_codes'] = kegg_enzyme_df['EC_substrates'].apply(lambda x: extract_compound_codes(','.join(x)) if isinstance(x, list) else [])
                kegg_enzyme_df['EC_product_codes'] = kegg_enzyme_df['EC_products'].apply(lambda x: extract_compound_codes(','.join(x)) if isinstance(x, list) else [])
                kegg_enzyme_df = kegg_enzyme_df.explode("EC_reactions")
                kegg_enzyme_df.to_pickle("kegg_enzyme_df.pkl")

            kegg_enzyme_df = kegg_enzyme_df.explode("EC_reactions")
            reactions = kegg_enzyme_df.EC_reactions.dropna().unique()
            progress.update(get_reactions, total = len(reactions))

            try:
                kegg_reaction_df = pd.read_pickle("kegg_reaction_df.pkl")
                progress.update(get_reactions, advance=len(reactions))
                time.sleep(0.1)
            except:
                kegg_reaction_dictionary = {}
                for i in range(0, len(reactions), n):
                    chunk = reactions[i:i + n]
                    reaction_dictionary = get_kegg_reactions(chunk, progress.console)
                    kegg_reaction_dictionary.update(reaction_dictionary)
                    progress.update(get_reactions, advance=len(reactions[i:i + n]))
                    time.sleep(1)


                #this is checking the number of equations returned in each reaction.
                lengths = [len(record['reaction_equation']) for record in kegg_reaction_dictionary.values()]
                counter = Counter(lengths)

                kegg_reaction_df = pd.DataFrame(kegg_reaction_dictionary).T
                kegg_reaction_df.to_pickle("kegg_reaction_df.pkl")

            kegg_reaction_enzyme_df = kegg_enzyme_df.merge(kegg_reaction_df, left_on = "EC_reactions", right_on = "reaction_id", how = "left")
            kegg_reaction_enzyme_df["reaction_substrate_codes"] = kegg_reaction_enzyme_df["reaction_substrate_codes"].apply(lambda d: d if isinstance(d, list) else [])
            kegg_reaction_enzyme_df["reaction_product_codes"] = kegg_reaction_enzyme_df["reaction_product_codes"].apply(lambda d: d if isinstance(d, list) else [])
            kegg_reaction_enzyme_df["EC_substrate_codes"] = kegg_reaction_enzyme_df["EC_substrate_codes"].apply(lambda d: d if isinstance(d, list) else [])
            kegg_reaction_enzyme_df["EC_product_codes"] = kegg_reaction_enzyme_df["EC_product_codes"].apply(lambda d: d if isinstance(d, list) else [])
            kegg_reaction_enzyme_df["entities"] = kegg_reaction_enzyme_df["EC_substrate_codes"] + kegg_reaction_enzyme_df["EC_product_codes"] + kegg_reaction_enzyme_df["reaction_substrate_codes"] + kegg_reaction_enzyme_df["reaction_product_codes"]
            
            kegg_reaction_enzyme_df.to_pickle("kegg_reaction_enzyme_df.pkl")

            compound_codes = kegg_reaction_enzyme_df.entities.explode().dropna().unique()
            
            progress.update(get_products_and_substrates, total=len(compound_codes))
            try:
                kegg_compounds_df = pd.read_pickle('kegg_compounds_df.pkl')
                progress.update(get_products_and_substrates, advance=len(compound_codes))
                time.sleep(0.1)
            except:
                kegg_compounds_dict = {}
                for compound in compound_codes:
                    compound_dict = get_kegg_compound(compound)
                    kegg_compounds_dict[compound] = compound_dict
                    progress.update(get_products_and_substrates, advance=1)
                    time.sleep(1)

                kegg_compounds_df = pd.DataFrame(kegg_compounds_dict).T

                kegg_compounds_df["kegg_compound_id"] = kegg_compounds_df.index

                kegg_compounds_df.loc[kegg_compounds_df.mol.isna(),["mol_source", "compound_id"]] = kegg_compounds_df.loc[kegg_compounds_df.mol.isna()].apply(lambda x: extract_secondary_id(x["kegg_compound_id"], x["dbxrefs"]), axis = 1, result_type = "expand").values
                kegg_compounds_df.to_pickle("kegg_compounds_df.pkl")
                
    try:
        kegg_compounds_df = pd.read_pickle("chebi_kegg_compounds_df.pkl")
    except:

        relations = pd.read_csv("chebi_relation.tsv", sep = "\t")
        structures = pd.read_csv("chebi_structures.csv")
        compounds = pd.read_csv("chebi_compounds.tsv", sep = "\t")

        chebi_subset = kegg_compounds_df.loc[(kegg_compounds_df.mol_source == "ChEBI") & (kegg_compounds_df.compound_id.isna() == False), ["compound_id"]].drop_duplicates("compound_id").copy()

        chebi_subset["chebi_result"] = chebi_subset["compound_id"].apply(lambda compound_id: find_chebi_structure(compound_id, compounds, relations, structures))
        chebi_subset = chebi_subset.explode("chebi_result")

        chebi_subset[["chebi_id", "chebi_mol", "chebi_relationship"]] = chebi_subset['chebi_result'].apply(pd.Series)
        chebi_subset = chebi_subset.dropna(subset=['chebi_mol'])
        
        kegg_compounds_df = kegg_compounds_df.merge(chebi_subset, on = "compound_id", how = "left")
        kegg_compounds_df["mol"] = kegg_compounds_df["mol"].combine_first(kegg_compounds_df["chebi_mol"])
        kegg_compounds_df.loc[kegg_compounds_df.mol.isna(),["mol_source", "compound_id"]] = kegg_compounds_df.loc[kegg_compounds_df.mol.isna()].apply(lambda x: extract_secondary_id(x["kegg_compound_id"], x["dbxrefs"], current_db = x["mol_source"]), axis = 1, result_type = "expand").values #update the secondary database identifier where chebi search still does not match a structure
        kegg_compounds_df.to_pickle("chebi_kegg_compounds_df.pkl")
        
    try:
        kegg_compounds_df = pd.read_pickle("pubchem_kegg_compounds_df.pkl")
    except:

        pubchem_keggs = kegg_compounds_df.loc[(kegg_compounds_df.mol.apply(lambda x: isinstance(x, Chem.Mol) == False)) & (kegg_compounds_df.mol_source == "PubChem")].compound_id.values.tolist()
        pubchem_inchi = pcp.get_properties("Inchi", pubchem_keggs, namespace=u'cid', as_dataframe=True).reset_index()
        pubchem_inchi.rename(columns = {"InChI": "mol"}, inplace = True)
        pubchem_inchi["CID"] = pubchem_inchi["CID"].astype("str")
        
        pubchem_inchi.loc[pubchem_inchi.mol.isna() == False, "mol"] = pubchem_inchi.loc[pubchem_inchi.mol.isna() == False].apply(lambda x: Chem.MolFromInchi(x["mol"]), axis = 1)
        kegg_compounds_df = kegg_compounds_df.merge(pubchem_inchi[["CID","mol"]], left_on = "compound_id", right_on = "CID", how = "left", suffixes=("", "_new"))
        kegg_compounds_df["mol"] = kegg_compounds_df["mol_new"].combine_first(kegg_compounds_df["mol"])
        kegg_compounds_df.drop(["mol_new", "CID"], axis = 1, inplace = True)
        kegg_compounds_df.to_pickle("pubchem_kegg_compounds_df.pkl")
        
    try:
        kegg_compounds_df = pd.read_pickle("glytoucan_kegg_compounds_df.pkl")
    except:
        glytoucan_ids_kegg = kegg_compounds_df.loc[(kegg_compounds_df.mol.isna()) & (kegg_compounds_df.mol_source == "GlyTouCan"), "compound_id"].values.tolist()
        if len(glytoucan_ids_kegg) > 0:
            glycosmos_glycans_wurcs = pd.read_csv("glycosmos_glycans_wurcs.csv")

            glytoucan_wurcs_kegg = glycosmos_glycans_wurcs.loc[glycosmos_glycans_wurcs["GlyTouCan ID"].isin(glytoucan_ids_kegg)].copy()
            glytoucan_wurcs_kegg["glycoct"] = glytoucan_wurcs_kegg.WURCS.apply(lambda x: get_glycoct_from_wurcs(x))
            glytoucan_wurcs_kegg["csdb_linear"] = glytoucan_wurcs_kegg.glycoct.apply(lambda x: get_csdb_from_glycoct(x))
            glytoucan_wurcs_kegg["mol"] = glytoucan_wurcs_kegg.apply(lambda x: get_smiles_from_csdb(x["csdb_linear"]), axis = 1)

            kegg_compounds_df = kegg_compounds_df.merge(glytoucan_wurcs_kegg[["GlyTouCan ID","mol"]], left_on = "compound_id", right_on = "GlyTouCan ID", how = "left", suffixes=("", "_new"))
            kegg_compounds_df["mol"] = kegg_compounds_df["mol_new"].combine_first(kegg_compounds_df["mol"])
            kegg_compounds_df.drop(["mol_new", "GlyTouCan ID"], axis = 1, inplace = True)
        kegg_compounds_df.to_pickle("glytoucan_kegg_compounds_df.pkl")
            
    try:
        kegg_compound_reaction_enzyme_df = pd.read_pickle("kegg_compound_reaction_enzyme_df.pkl")
    except:
        kegg_reaction_enzyme_df = kegg_reaction_enzyme_df.explode("entities")
        kegg_compound_reaction_enzyme_df = kegg_reaction_enzyme_df.merge(kegg_compounds_df, left_on = "entities", right_on = "kegg_compound_id", how = "left")
    
        kegg_compound_reaction_enzyme_df.to_pickle("kegg_compound_reaction_enzyme_df.pkl")
    
    missing_entities = kegg_compound_reaction_enzyme_df.loc[((kegg_compound_reaction_enzyme_df.mol.isna()) & (kegg_compound_reaction_enzyme_df.EC_reactions.isna()) & 
                                (kegg_compound_reaction_enzyme_df.EC_substrate_codes.str.len() == 0) & 
                                (kegg_compound_reaction_enzyme_df.EC_product_codes.str.len() == 0) & 
                                (kegg_compound_reaction_enzyme_df.matched_name != "Deleted entry")) | ((kegg_compound_reaction_enzyme_df.mol.isna()) & 
        (kegg_compound_reaction_enzyme_df.matched_name != "Deleted entry")),  ["matched_entry","reaction_text"]].drop_duplicates(["matched_entry"]).dropna(subset = ["reaction_text"]).copy()
    
    

    try:
        missing_entities_solved = pd.read_pickle("missing_entities_final.pkl")
    except:
        os.environ['OPENAI_API_KEY'] = args.openai_api_key
        os.environ['OPENAI_API_ORG'] = args.openai_api_org
        
        nlp = assemble("spacy_config.cfg")

        relations = pd.read_csv("chebi_relation.tsv", sep = "\t")
        structures = pd.read_csv("chebi_structures.csv")
        compounds = pd.read_csv("chebi_compounds.tsv", sep = "\t")

        missing_entities = missing_entities.explode("reaction_text").drop_duplicates(subset = ["reaction_text"]).reset_index(drop = True)
        missing_entities = missing_entities.dropna(subset = "reaction_text")
        missing_entities["reaction_text"] = missing_entities["reaction_text"].str.strip()
        missing_entities["reaction_text"] = missing_entities.reaction_text.str.replace('^\(\d+[A-Za-z]*\) ', '', regex = True)
        missing_entities["reaction_text"] = missing_entities.reaction_text.str.replace(' *an ', '', regex = True)
        missing_entities[["lhs", "rhs"]] = missing_entities.reaction_text.str.split("=", n=1, expand=True)
        missing_entities["lhs"] = missing_entities.lhs.str.split("+")
        missing_entities["rhs"] = missing_entities.rhs.str.split("+")

        missing_entities = missing_entities.explode("lhs")
        missing_entities = missing_entities.explode("rhs")
        missing_entities.reset_index(drop = True, inplace = True)

        missing_entities["lhs"] = missing_entities.lhs.str.replace(" *\d+ ", "", regex = True).str.strip()
        missing_entities["rhs"] = missing_entities.rhs.str.replace(" *\d+ ", "", regex = True).str.strip()
        missing_entities[["lhs", "rhs"]] = missing_entities[["lhs", "rhs"]].replace(' ', np.nan)

        missing_entities.loc[missing_entities.rhs.isna(), "search_type"] = "NER"
        missing_entities.loc[missing_entities.search_type.isna(), "search_type"] = "entities"

        missing_entities.loc[missing_entities.search_type == "entities", "entities"] = missing_entities.loc[missing_entities.search_type == "entities"].apply(lambda row: [value for value in [row['lhs'], row['rhs']] if value is not None and not pd.isna(value)], axis=1)

        missing_entities.loc[missing_entities.search_type == "NER", "entities"] = missing_entities.loc[missing_entities.search_type == "NER"].apply(lambda x: get_spacy_entities(nlp, x["reaction_text"]), axis = 1)

        missing_entities.to_pickle("missing_entities.pkl")

        missing_entities_entities_unique = missing_entities[["entities"]].drop_duplicates(subset = "entities").dropna()

        missing_entities_entities_unique["chebi_result"] = missing_entities_entities_unique.entities.apply(lambda x: get_compound_id_by_name(compounds, x, compounds, relations, structures, search = True))

        missing_entities_entities_unique = missing_entities_entities_unique.explode("chebi_result")

        missing_entities_entities_unique[["chebi_id", "chebi_mol", "chebi_relationship"]] = missing_entities_entities_unique["chebi_result"].apply(lambda x: pd.Series(x) if isinstance(x, list) else pd.Series([np.nan, np.nan, np.nan]))

        missing_entities_entities_unique = missing_entities_entities_unique.loc[(missing_entities_entities_unique.chebi_mol.isna() == False) & (missing_entities_entities_unique.chebi_mol != "no structure identified")]

        missing_entities_solved = missing_entities.merge(missing_entities_entities_unique, on = "entities", how = "left")

        missing_entities_solved = missing_entities_solved.loc[missing_entities_solved.chebi_mol.isna() == False]

        missing_entities_solved.to_pickle("missing_entities_final.pkl")

    kegg_compound_reaction_enzyme_df = kegg_compound_reaction_enzyme_df.merge(missing_entities_solved, on = "matched_entry", how = "left", suffixes=("", "_new"))

    kegg_compound_reaction_enzyme_df["mol"] = kegg_compound_reaction_enzyme_df["chebi_mol_new"].combine_first(kegg_compound_reaction_enzyme_df["mol"])
    kegg_compound_reaction_enzyme_df["chebi_result"] = kegg_compound_reaction_enzyme_df["chebi_result_new"].combine_first(kegg_compound_reaction_enzyme_df["chebi_result"])
    kegg_compound_reaction_enzyme_df["chebi_id"] = kegg_compound_reaction_enzyme_df["chebi_id_new"].combine_first(kegg_compound_reaction_enzyme_df["chebi_id"])
    kegg_compound_reaction_enzyme_df["chebi_relationship"] = kegg_compound_reaction_enzyme_df["chebi_relationship_new"].combine_first(kegg_compound_reaction_enzyme_df["chebi_relationship"])
    kegg_compound_reaction_enzyme_df.drop(columns = ["chebi_mol_new", "chebi_result_new", "chebi_id_new", "chebi_relationship_new", "reaction_text"], inplace = True)

    kegg_compound_reaction_enzyme_df['compound_id'].fillna(kegg_compound_reaction_enzyme_df['chebi_id'], inplace=True)
    kegg_compound_reaction_enzyme_df["unique_id"] = kegg_compound_reaction_enzyme_df["mol_source"] + "_" + kegg_compound_reaction_enzyme_df["compound_id"].astype("str")
    kegg_compound_reaction_enzyme_df.loc[kegg_compound_reaction_enzyme_df.chebi_id.isna() == False, "unique_id"] = "Chebi" + "_" + kegg_compound_reaction_enzyme_df["chebi_id"].astype("str")
    kegg_compound_reaction_enzyme_df.to_pickle("final_kegg_compound_reaction_enzyme_df.pkl")
