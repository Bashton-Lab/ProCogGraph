#!/usr/bin/env python

import argparse
import pandas as pd
from gemmi import cif
from utils import process_ec_records, get_updated_enzyme_records, get_scop_domains_info, extract_interpro_domain_annotations, get_pfam_annotations, parse_cddf, build_cath_dataframe
import numpy as np
from Bio.ExPASy import Enzyme as EEnzyme
import re

#copied from extract_pdbe_info for now - combine these into a file that they can be used from for snakemake and nextflow pipelines
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
import requests
import json
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

def get_sugar_smiles_from_wurcs(wurcs_list, csdb_linear_cache, smiles_cache, glycoct_cache):
    #work out how to update the caches here - needs to return some updated caches from the func to save in main
    sugar_smiles = {}
    for wurcs in wurcs_list:
        smiles = None
        glycoct = get_glycoct_from_wurcs(wurcs, glycoct_cache)
        if glycoct is not None:
            csdb = get_csdb_from_glycoct(glycoct, csdb_linear_cache)
            if csdb is not None:
                smiles = get_smiles_from_csdb(csdb, smiles_cache)
        sugar_smiles[wurcs] = smiles
    return sugar_smiles

def get_chem_comp_descriptors(ccd_doc, comp_id_list):
    ligand_descriptors = {}
    for ligand in comp_id_list:
        lig_descriptor = None
        lig_block = ccd_doc.find_block(ligand)
        lig_descriptors = pd.DataFrame(lig_block.find_mmcif_category("_pdbx_chem_comp_descriptor."), columns = ["comp_id", "type", "program", "program_version", "descriptor"])
        lig_descriptors["descriptor"] = lig_descriptors.descriptor.str.strip("\"|'")
        lig_descriptor = lig_descriptors.loc[lig_descriptors.type == "SMILES"].descriptor.values[0]
        ligand_descriptors[ligand] = lig_descriptor
    return ligand_descriptors

def process_sifts_ec_map(sifts_ec_mapping_file, ec_records_file):
    sifts_chains_ec = sifts_ec_mapping_file.loc[sifts_ec_mapping_file.EC_NUMBER != "?"].copy()
    sifts_chains_uniprot = sifts_ec_mapping_file.groupby(["PDB", "CHAIN"]).agg({"ACCESSION": set}).reset_index().copy()
    sifts_chains_uniprot["ACCESSION"] = sifts_chains_uniprot["ACCESSION"].apply(lambda x: "|".join(x)) #join the list of uniprot accessions with a pipe for downstream neo4j integration
    sifts_chains_uniprot.rename(columns = {"ACCESSION" : "uniprot_accession"}, inplace = True)

    sifts_chains_ec = sifts_chains_ec.groupby(["PDB"]).agg({"EC_NUMBER": set}).reset_index() #group these into a set of pdb associated ec's
    sifts_chains_ec["EC_NUMBER"] = sifts_chains_ec["EC_NUMBER"].apply(lambda x: ",".join(x)) #join the list of EC numbers into a single string for ec records function 
    sifts_chains_ec = get_updated_enzyme_records(sifts_chains_ec, ec_records_file, ec_col = "EC_NUMBER")
    sifts_chains_ec.rename(columns = {"EC_NUMBER": "protein_entity_ec", "ACCESSION" : "uniprot_accession"}, inplace = True)

    return sifts_chains_ec, sifts_chains_uniprot

def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--ccd_cif', type=str, help='cif file containing chemical component dictionary in mmcif format')
    parser.add_argument('--contacts_file', type=str, help='csv file containing concatenated contacts from pdbe-arpeggio for all pdbs')
    parser.add_argument('--pfam_a_file', type=str, help='pfam a file')
    parser.add_argument('--pfam_clan_rels', type=str, help='pfam clan relationships file')
    parser.add_argument('--pfam_clans', type=str, help='pfam clans file')
    parser.add_argument('--scop_domains_info_file', type=str, help='scop domains info file')
    parser.add_argument('--scop_descriptions_file', type=str, help='scop descriptions file')
    parser.add_argument('--interpro_xml', type=str, help='interpro xml file')
    parser.add_argument('--cddf', type=str, help='cath domain description file')
    parser.add_argument('--glycoct_cache', type=str, help='glycoct cache file')
    parser.add_argument('--smiles_cache', type=str, help='smiles cache file')
    parser.add_argument('--csdb_linear_cache', type=str, help='csdb linear cache file')
    parser.add_argument('--enzyme_dat_file', type=str, help='enzyme.dat file')
    parser.add_argument('--enzyme_class_file', type=str, help='enzclass.txt file')
    parser.add_argument('--sifts_ec_mapping', type=str, help='sifts ec mapping file')

    args = parser.parse_args()  
    
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

    contacts = pd.read_csv(args.contacts_file, sep = "\t", na_values = ["NaN", "None"], keep_default_na = False)   #for occurences of sodium

    ccd_doc = cif.read(args.ccd_cif)
    
    ligand_ids = contacts.loc[contacts.type == "ligand"].hetCode.unique()
    ligand_descriptors = get_chem_comp_descriptors(ccd_doc, ligand_ids)
    contacts.loc[contacts.type == "ligand", "descriptor"] = contacts.loc[contacts.type == "ligand"]["hetCode"].apply(lambda x: ligand_descriptors[x])

    glycoct_cache = pd.read_pickle(f"{args.glycoct_cache}")
    smiles_cache = pd.read_pickle(f"{args.smiles_cache}")   
    csdb_linear_cache = pd.read_pickle(f"{args.csdb_linear_cache}")

    wurcs_list = contacts.loc[contacts.type == "sugar", "descriptor"].unique()
    sugar_smiles = get_sugar_smiles_from_wurcs(wurcs_list, csdb_linear_cache, smiles_cache, glycoct_cache)

    contacts.loc[contacts.type == "sugar", "descriptor"] = contacts.loc[contacts.type == "sugar"]["descriptor"].apply(lambda x: sugar_smiles[x])
    #after merging the contacts with descriptors we probably want to flag or remove the ones with failed descriptors
    ##then we need to assign a ligand uniqueID to all unique ligands in the contacts file - check implementation of this
    
    #assign ec and uniprot information to contacts

    ec_records_df_grouped = process_ec_records(args.enzyme_dat_file, args.enzyme_class_file)
    ec_records_df = ec_records_df_grouped.explode("ID")

    sifts_chains = pd.read_csv(f"{args.sifts_ec_mapping}", sep = "\t", comment="#")
    sifts_chains_ec, sifts_chains_uniprot = process_sifts_ec_map(sifts_chains, ec_records_df)
    
    contacts_ec = contacts.merge(sifts_chains_ec, left_on = "pdb_id", right_on = "PDB", how = "left", indicator = True)
    contacts_ec_unmatched = contacts_ec.loc[contacts_ec._merge != "both"].copy().drop(columns = ["_merge", "PDB"])
    contacts_ec = contacts_ec.loc[contacts_ec._merge == "both"].copy().drop(columns = ["_merge", "PDB"])
    contacts_ec_unmatched.to_csv(f"contacts_ec_unmatched.csv", index = False)
    
    contacts_ec = contacts_ec.merge(sifts_chains_uniprot, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left")
    contacts_ec.drop(columns = ["PDB", "CHAIN"], inplace = True)
    #contacts_ec.to_csv(f"combined_contacts_processed.tsv", sep = "\t" , index = False)
    #now add uniprot info to the contacts file
    contacts_ec_uniprot = contacts_ec.merge(sifts_chains_uniprot, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left", indicator = True)
    #what are we validating here to need the indicator col
    contacts_ec_uniprot.drop(columns = ["PDB", "CHAIN", "_merge"], inplace = True)
    
    #get the unique ligands from the dataset, and assign uniqueIDs - extract unique ligands to file for scoring.

    bound_entities_to_score = contacts_ec_uniprot[["descriptor", "description", "hetCode", "type", "ec_list"]].groupby(["hetCode", "description", "descriptor"]).agg({"ec_list": set}).reset_index().reset_index()
    print(bound_entities_to_score)
    #add ligand entity id to the contacts file
    contacts_ec_uniprot = contacts_ec_uniprot.merge(bound_entities_to_score[["descriptor", "description", "hetCode", "index"]], on = ["descriptor", "description", "hetCode"], how = "left", indicator = True)
    assert(len(contacts_ec_uniprot.loc[contacts_ec_uniprot._merge != "both"]) == 0)
    contacts_ec_uniprot.drop(columns = "_merge", inplace = True)
    contacts_ec_uniprot.rename(columns = {"index": "ligand_uniqueID"}, inplace = True)

    bound_entities_to_score.rename(columns = {"index" : "ligand_entity_id", "hetCode": "bl_name"}, inplace = True)
    bound_entities_to_score.to_pickle(f"bound_entities_to_score.pkl")

    

    #split the domain information into separate database specific dataframes for annotation
    cath_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "CATH"].copy()
    scop_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "SCOP"].copy()
    scop2b_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "SCOP2B"].copy()
    pfam_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "Pfam"].copy()
    interpro_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "InterPro"].copy()

    core_cols = contacts_ec_uniprot.columns.tolist()
    cath_cols = ['cath_db_version', 'cath_db_verdate', 'cath_name', 'cath_source', 'cath_code', 'cath_class', 'cath_architecture', 'cath_topology',
        'cath_homologous_superfamily', 'cath_domain_length', 'cath_domain_seq_header', 'cath_domain_seqs', 'cath_num_segments', 'cath_segments_dict']
    scop_cols = ['scop_id', 'sccs', 'domain_sunid', 'ancestor_sunid', 'cl_id', 'cf_id', 'sf_id', 'fa_id', 'dm_id', 'sp_id',
        'px_id', 'cl_description', 'cf_description', 'sf_description','fa_description', 'dm_description', 'sp_description', 'px_description']
    scop2b_cols = []
    pfam_cols = ['clan', 'clan_acc', 'clan_comment', 'clan_description', 'pfam', 'pfam_accession', 'pfam_description', 'pfam_name']
    interpro_cols = ['interpro_id','interpro_short_name', 'dbxref']


    #process domain database information
    if len(cath_contacts) > 0:
        #run the cath_domain_list on only the domains we need
        cath_domain_list = cath_contacts.xref_db_acc.unique()
        cath_parsed_data = parse_cddf(args.cddf, cath_domain_list)
        cath_domains_info = build_cath_dataframe(cath_parsed_data)
        cath_contacts = cath_contacts.merge(cath_domains_info, how = "left", left_on = "xref_db_acc", right_on = "cath_domain", indicator = True)
        print(cath_contacts.loc[cath_contacts._merge != "both", ["xref_db_acc", "cath_domain"]])
        assert(len(cath_contacts.loc[cath_contacts._merge != "both"]) == 0)
        cath_contacts.drop(columns = "_merge", inplace = True)
        cath_contacts.to_csv(f"cath_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        cath_contacts = pd.DataFrame(columns = core_cols + cath_cols)
        cath_contacts.to_csv(f"cath_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")

    if len(scop_contacts) > 0:
        scop_domains_info = get_scop_domains_info(args.scop_domains_info_file, args.scop_descriptions_file)
        scop_contacts = scop_contacts.merge(scop_domains_info, how = "left", left_on = "xref_db_acc", right_on = "scop_id", indicator = True)
        assert(len(scop_contacts.loc[scop_contacts._merge != "both"]) == 0)
        scop_contacts.drop(columns = "_merge", inplace = True)
        scop_contacts.to_csv(f"scop_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        scop_contacts = pd.DataFrame(columns = core_cols + scop_cols)
        scop_contacts.to_csv(f"scop_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")

    if len(pfam_contacts) > 0:
        pfam_annotations = get_pfam_annotations(args.pfam_a_file, args.pfam_clan_rels, args.pfam_clans)
        pfam_contacts = pfam_contacts.merge(pfam_annotations, left_on = "xref_db_acc" , right_on = "pfam_accession", how = "left")
        pfam_contacts.to_csv(f"pfam_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        pfam_contacts = pd.DataFrame(columns = core_cols + pfam_cols)
        pfam_contacts.to_csv(f"pfam_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")

    if len(interpro_contacts) > 0:
        interpro_annotations = extract_interpro_domain_annotations(args.interpro_xml)
        interpro_contacts = interpro_contacts.merge(interpro_annotations, left_on = "xref_db_acc", right_on = "interpro_accession", how = "left")
        interpro_contacts = interpro_contacts.loc[(interpro_contacts.xref_db_acc.isna() == False) & ((interpro_contacts.xref_db_acc.str.contains("SSF")) | (interpro_contacts.xref_db_acc.str.contains("G3DSA")))]
        superfamily_contacts = interpro_contacts.loc[(interpro_contacts.xref_db_acc.isna() == False) & (interpro_contacts.xref_db_acc.str.contains("SSF"))]
        g3dsa_contacts = interpro_contacts.loc[(interpro_contacts.xref_db_acc.isna() == False) & (interpro_contacts.xref_db_acc.str.contains("G3DSA"))]
        interpro_contacts.to_csv(f"interpro_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
        superfamily_contacts.to_csv(f"superfamily_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
        gene3d_sa_contacts.to_csv(f"g3dsa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        interpro_contacts = pd.DataFrame(columns = core_cols + interpro_cols)
        superfamily_contacts = pd.DataFrame(columns = core_cols + interpro_cols)
        gene3d_sa_contacts = pd.DataFrame(columns = core_cols + interpro_cols)
        interpro_contacts.to_csv(f"interpro_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
        superfamily_contacts.to_csv(f"superfamily_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
        gene3d_sa_contacts.to_csv(f"g3dsa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    #output the processed contacts files
    if len(scop2b_contacts) > 0:
        scop2b_contacts.to_csv(f"scop2b_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        scop2b_contacts = pd.DataFrame()
        scop2b_contacts.to_csv(f"scop2b_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")

if __name__ == "__main__":
    main()