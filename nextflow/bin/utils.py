#!/usr/bin/env python

import re
import numpy as np
import pandas as pd
import requests
from urllib.parse import quote
from bs4 import BeautifulSoup
import json 
from Bio.ExPASy import Enzyme as EEnzyme
from pdbeccdutils.helpers.mol_tools import fix_molecule
from rdkit import Chem
import gzip 
import xml.etree.ElementTree as ET

#make this function be applicable to extract_pdbe_info script too - need to check if the grouped output at end is appropriate.
def process_ec_records(enzyme_dat_file, enzyme_class_file):
    """
    Process EC records and generate related dataframes.

    Returns:
    ec_records_df_grouped (DataFrame): Grouped EC records dataframe.
    ec_class_descriptions (DataFrame): EC class descriptions dataframe.
    ec_subclass_descriptions (DataFrame): EC subclass descriptions dataframe.
    ec_subsubclass_descriptions (DataFrame): EC subsubclass descriptions dataframe.
    """

    with open(enzyme_dat_file) as handle:
        ec_records = EEnzyme.parse(handle)
        ec_records_list = []
        for record in ec_records: 
            ec_record_series = pd.Series(record)
            ec_records_list.append(ec_record_series)

    ec_records_df = pd.DataFrame(ec_records_list)
    ec_records_df["TRANSFER"] = ec_records_df.apply(lambda x: get_terminal_record(x["ID"], x, ec_records_df), axis = 1)
    ec_records_df["TRANSFER"] = ec_records_df["TRANSFER"].fillna(ec_records_df.ID)

    with open(enzyme_class_file, 'r') as file:
        lines = [line.strip() for line in file if re.match(r'^\d\.', line.strip())]

    ec_descriptions = []
    for line in lines:
        ec = ''.join(line[0:9]).replace(' ', '')
        description = re.findall(r'.*.-\s+(.*)', line)[0]
        ec_descriptions.append((ec, description))

    ec_descriptions_df = pd.DataFrame(ec_descriptions, columns= ["EC", "description"])
    ec_records_df_grouped = ec_records_df.groupby("TRANSFER").agg({"ID": list}).reset_index()

    ec_records_df_grouped = ec_records_df_grouped.merge(ec_records_df[["ID", "DE"]], left_on = "TRANSFER", right_on = "ID")
    ec_records_df_grouped.rename(columns = {"ID_x" : "ID"}, inplace = True)
    ec_records_df_grouped.drop(columns = ["ID_y"], inplace = True)
    ec_records_df_grouped[["class", "subclass", "subsubclass", "term"]] = ec_records_df_grouped.TRANSFER.str.split(".", expand = True)
    ec_records_df_grouped = ec_records_df_grouped.loc[ec_records_df_grouped.DE != "Deleted entry."].copy()
    
    subsubclass = ec_records_df_grouped["class"].astype("str") + "." + ec_records_df_grouped["subclass"] + "." + ec_records_df_grouped["subsubclass"] + ".-"
    subclass = ec_records_df_grouped["class"].astype("str") + "." + ec_records_df_grouped["subclass"] + ".-.-"
    ec_class = ec_records_df_grouped["class"].astype("str") + ".-.-.-"

    ec_records_df_grouped["subsubclass"] = subsubclass
    ec_records_df_grouped["subclass"] = subclass
    ec_records_df_grouped["class"] = ec_class

    ec_class_descriptions = ec_descriptions_df.loc[ec_descriptions_df.EC.str.contains("\d\.-\.-\.-$", regex = True)].copy().rename(columns = {"EC": "class","description" : "class_description"})
    ec_subclass_descriptions = ec_descriptions_df.loc[ec_descriptions_df.EC.str.contains(r"\d\.-\.-$", regex = True)].copy().rename(columns = {"EC":"subclass","description" : "subclass_description"})
    ec_subsubclass_descriptions = ec_descriptions_df.loc[ec_descriptions_df.EC.str.contains(r"\d\.-$", regex = True)].copy().rename(columns = {"EC" : "subsubclass","description" : "subsubclass_description"})

    ec_records_df_grouped = ec_records_df_grouped.merge(ec_class_descriptions, on = "class", how = "left")
    ec_records_df_grouped = ec_records_df_grouped.merge(ec_subclass_descriptions, on = "subclass", how = "left")
    ec_records_df_grouped = ec_records_df_grouped.merge(ec_subsubclass_descriptions, on = "subsubclass", how = "left")

    return ec_records_df_grouped

#this is a core function - should be stored elsewhere.
#a limitation of this function is that it selects only one representative EC number for transferred entries
def get_terminal_record(entry, row, df):
    entry = row.ID
    pattern = r"[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+"
    while row.DE.startswith("Transferred entry: "):
        transfers = re.findall(pattern, row.DE)
        #when updating, if multiple possible transfers, selecting only first.
        row = df.loc[df.ID == transfers[0]].iloc[0]
    return row.ID

def get_csdb_from_glycoct(glycoct, cache_df):
    if glycoct is np.nan or glycoct == None:
        return np.nan
    elif glycoct in cache_df.glycoct.values:
        csdb_linear = cache_df.loc[cache_df.glycoct == glycoct, "csdb"].values[0]
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
    
def get_glycoct_from_wurcs(wurcs, cache_df):
    if wurcs is np.nan or wurcs == None:
        return np.nan
    elif wurcs in cache_df.WURCS.values:
        glycoct = cache_df.loc[cache_df.WURCS == wurcs, "glycoct"].values[0]
        return glycoct
    else:
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

def get_smiles_from_csdb(csdb_linear, cache_df):
    if csdb_linear is np.nan or csdb_linear == None:
        return np.nan
    elif csdb_linear in cache_df.csdb.values:
        smiles = cache_df.loc[cache_df.csdb == csdb_linear, "descriptor"].values[0]
        return smiles
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
    
def pdbe_sanitise_smiles(smiles, return_mol = False, return_sanitisation = False):
    """
    Sanitises a smiles string using pdbeccdutils fix_molecule functions and
    returns a canonical smiles string or RDKit molecule object. Requires that
    smiles string can be loaded into an RDKit molecule object - returns none if
    this is not possible.
    """
    try:
        mol = Chem.MolFromSmiles(Chem.CanonSmiles(smiles))
    except:
        if return_sanitisation:
            return None, None
        else:
            return None
    
    sanitisation = fix_molecule(mol)
    if sanitisation:
        if return_mol:
            if return_sanitisation:
                return mol, sanitisation
            else:
                return mol
        else:
            sanitised_smiles = Chem.CanonSmiles(Chem.MolToSmiles(mol))
            if return_sanitisation:
                return sanitised_smiles, sanitisation
            else:
                return sanitised_smiles
    else:
        sanitised_smiles = None
        if return_sanitisation:
            return sanitised_smiles, sanitisation
        else:
            return sanitised_smiles
        
def extract_interpro_domain_annotations(xml_file):
    with gzip.open(xml_file, 'rb') as f:
        tree = ET.parse(f)
        root = tree.getroot()

        # Find all interpro elements in the XML file
        interpro_elements = root.findall('.//interpro')

        interpro_info = []
        # Iterate through each interpro element
        for interpro in interpro_elements:
            superfamily_annotations = {}
            gene3d_annotations = {}
            interpro_id = interpro.attrib['id']
            interpro_short_name = interpro.attrib['short_name']
            # Store PFAM annotations for the interpro ID
            interpro_info.append({"interpro_accession": interpro_id,
                                        "interpro_name": interpro_short_name})
    interpro_info_df = pd.DataFrame(interpro_info)       

    return interpro_info_df

def get_scop_domains_info(domain_info_file, descriptions_file):
    
    def clean_and_merge_scop_col(df, column_id, description_df):
        level = df[column_id].str.split("=").str.get(0).values[0]
        df[column_id] = df[column_id].str.split("=").str.get(1).astype(int)
        df = df.merge(description_df.loc[description_df.level == level, ["level_sunid", "level", "level_description"]],left_on = column_id, right_on = "level_sunid", indicator = True)
        df.rename(columns = {"level_description": f"{level}_description"}, inplace = True)
        assert len(df.loc[df._merge != "both"]) == 0
        df.drop(columns = ["_merge", "level_sunid", "level"], inplace = True)
        return df
    
    scop_domains_info = pd.read_csv(domain_info_file, sep = "\t", comment = "#", header = None, names = ["scop_id", "pdb_id", "scop_description", "sccs", "domain_sunid", "ancestor_sunid"])
    scop_id_levels = ["cl_id", "cf_id", "sf_id", "fa_id", "dm_id", "sp_id", "px_id"]
    scop_domains_info[scop_id_levels] = scop_domains_info.ancestor_sunid.str.split(",", expand = True)
    scop_descriptions = pd.read_csv(descriptions_file, sep = "\t", comment = "#" , header = None, names = ["level_sunid", "level", "level_sccs", "level_sid", "level_description"])

    for column in scop_id_levels:
        scop_domains_info = clean_and_merge_scop_col(scop_domains_info, column, scop_descriptions)
    
    scop_domains_info.drop(columns = ["pdb_id", "scop_description"], inplace = True)
    return scop_domains_info

def get_pfam_annotations(pfam_a_file, clan_membership_file, clan_info_file):
    pfam_a = pd.read_csv(pfam_a_file, sep = "\t", header = None, usecols = [0,1,3], names = ["pfam_accession", "pfam_name", "pfam_description"])
    pfam_clan_rels = pd.read_csv(clan_membership_file, sep = "\t", header = None, names = ["clan", "pfam"])
    pfam_clans = pd.read_csv(clan_info_file, sep = "\t", comment = "#", header = None, names = ["clan_acc", "clan_id", "previous_id", "clan_description", "clan_author", "deposited_by", "clan_comment", "updated", "created", "version", "number_structures", "number_archs", "number_species", "number_sequences", "competed", "uniprot_competed"])
    pfam_clan_df = pfam_clan_rels.merge(pfam_clans[["clan_acc", "clan_description", "clan_comment"]], left_on = "clan", right_on = "clan_acc", how = "left", indicator = True)
    assert(len(pfam_clan_df.loc[pfam_clan_df._merge != "both"]) == 0)
    pfam_clan_df.drop(columns = "_merge", inplace = True)
    
    pfam_a_clans_merged = pfam_a.merge(pfam_clan_df, left_on = "pfam_accession", right_on = "pfam", how = "left")
    return pfam_a_clans_merged

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
    residue_ec_records_grouped["TRANSFER"] = residue_ec_records_grouped["TRANSFER"].str.strip(",")
    residue_ec_records_grouped.rename(columns = {"TRANSFER" : "ec_list"}, inplace = True)
    
    df_merged = df.merge(residue_ec_records_grouped, on = ec_col, how = "left", indicator = True)
    assert(len(df_merged.loc[df_merged["_merge"] != "both"]) == 0)
    df_merged.drop(columns = "_merge", inplace = True)
    df_merged = df_merged.loc[df_merged.ec_list != ""] #remove any rows where the ec_list is empty - we cant process these anyway.
    return(df_merged)

def parse_cddf(file_path, domain_list):
    with open(file_path, 'r') as file:
        entries = {}
        lines = []
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                continue
            if line.startswith("//"):
                entries[domain] = lines
                lines = []
            elif line.startswith("DOMAIN"):
                domain = line[10:]
                lines.append(line)
            else:
                lines.append(line)
    matching_entries = []
    for domain in domain_list:
        if domain in entries.keys():
            current_entry = {}
            current_segment = {}
            for line in entries[domain]:
                if line.startswith('SEGMENT'):
                    if 'SEGMENTS' not in current_entry:
                        current_entry['SEGMENTS'] = []
                    if current_segment:
                        current_entry['SEGMENTS'].append(current_segment)
                    current_segment = {'SEGMENT': line.split()[1]}
                elif line.startswith('ENDSEG'):
                    if current_segment:
                        current_entry['SEGMENTS'].append(current_segment)
                    current_segment = {}
                else:
                    tag = line[0:9].strip()
                    value = line[10:]

                    if tag == 'FORMAT':
                        if 'FORMAT' not in current_entry:
                            current_entry['FORMAT'] = []
                        current_entry['FORMAT'].append(value)
                    elif tag in ['DOMAIN', 'PDBID', 'VERSION', 'VERDATE', 'NAME', 'SOURCE', 'CATHCODE', 'CLASS', 'ARCH', 'TOPOL', 'HOMOL', 'DLENGTH', 'DSEQH', 'DSEQS', 'NSEGMENTS']:
                        if tag not in current_entry:
                            current_entry[tag] = []
                        current_entry[tag].append(value)
                    elif tag in ['SRANGE', 'SLENGTH', 'SSEQH', 'SSEQS']:
                        current_segment[tag] = value
                    elif tag == 'COMMENTS':
                        if 'COMMENTS' not in current_entry:
                            current_entry['COMMENTS'] = []
                        current_entry['COMMENTS'].append(value)
                    else:
                        # Invalid tag
                        continue
            matching_entries.append(current_entry)
    return matching_entries


def build_cath_dataframe(parsed_data):
    combine_keys = ["DOMAIN", "PDBID", "VERSION", "VERDATE", "NAME", "SOURCE", "CATHCODE", "CLASS", "ARCH", "TOPOL", "HOMOL", "DLENGTH", "DSEQH", "DSEQS", "NSEGMENTS"]
    cath_df_columns = {"DOMAIN": "cath_domain" , "VERSION": "cath_db_version", "VERDATE": "cath_db_verdate", "NAME":"cath_name", "SOURCE": "cath_source", "CATHCODE": "cath_code", "CLASS": "cath_class_name", "ARCH":"cath_architecture_name", "TOPOL": "cath_topology_name", "HOMOL": "cath_homologous_superfamily_name", "DLENGTH": "cath_domain_length", "DSEQH":"cath_domain_seq_header", "DSEQS": "cath_domain_seqs", "NSEGMENTS": "cath_num_segments", "SEGMENTS": "cath_segments_dict"}
    dfs = []
    for entry in parsed_data:
        df_dict = {}
        for key, values in entry.items():
            if key in combine_keys:
                df_dict[cath_df_columns[key]] = ' '.join(values) if isinstance(values, list) else values
            elif key == "FORMAT":
                continue
            else:
                df_dict[cath_df_columns[key]] = values
        dfs.append(pd.DataFrame([df_dict]))
    dfs_combined = pd.concat(dfs, ignore_index=True)
    topology_regex = r'^(\d+\.\d+\.\d+)\.'
    architecture_regex = r'^(\d+\.\d+)\.'
    class_regex = r'^(\d+)\.'
    
    dfs_combined["cath_class"] = dfs_combined.cath_code.str.extract(class_regex, expand = True)
    dfs_combined["cath_architecture"] = dfs_combined.cath_code.str.extract(architecture_regex, expand = True)
    dfs_combined["cath_topology"] = dfs_combined.cath_code.str.extract(topology_regex, expand = True)
    dfs_combined["cath_homologous_superfamily"] = dfs_combined["cath_code"]
    return dfs_combined

def build_g3dsa_dataframe(cath_names, cath_domain_list): #unlike cath, g3dsa annotations are at superfamily level (also works for xml cath annotations from sifts which are at homolsuperfam level)
    topology_regex = r'^(\d+\.\d+\.\d+)\.'
    architecture_regex = r'^(\d+\.\d+)\.'
    class_regex = r'^(\d+)\.'
    domain_list = []
    for homologous_superfamily in cath_domain_list:
        homologous_superfamily_name = cath_names.loc[cath_names.cath_code == homologous_superfamily, "name"].values[0]
        topology = re.search(topology_regex, homologous_superfamily).group(1)
        topology_name = cath_names.loc[cath_names.cath_code == topology, "name"].values[0]
        architecture = re.search(architecture_regex, homologous_superfamily).group(1)
        architecture_name = cath_names.loc[cath_names.cath_code == architecture, "name"].values[0]
        class_ = re.search(class_regex, homologous_superfamily).group(1)
        class_name = cath_names.loc[cath_names.cath_code == class_, "name"].values[0]
        domain = {"cath_domain": homologous_superfamily, "cath_name" : homologous_superfamily_name, "cath_code": homologous_superfamily, "cath_homologous_superfamily": homologous_superfamily, "cath_homologous_superfamily_name": homologous_superfamily_name, "cath_topology": topology, "cath_topology_name": topology_name,  "cath_architecture": architecture, "cath_architecture_name": architecture_name, "cath_class": class_, "cath_class_name" : class_name}
        domain_list.append(domain)
    domain_df = pd.DataFrame(domain_list)
    return domain_df