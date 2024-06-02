#!/usr/bin/env python

import argparse
import pandas as pd
from gemmi import cif
from utils import process_ec_records, get_updated_enzyme_records, get_scop_domains_info, extract_interpro_domain_annotations, get_pfam_annotations, get_glycoct_from_wurcs, get_csdb_from_glycoct, get_smiles_from_csdb, build_cath_dataframe, parse_cddf, build_g3dsa_dataframe
import numpy as np
from Bio.ExPASy import Enzyme as EEnzyme
import re
from urllib.parse import quote
from rdkit.Chem import PandasTools

def get_sugar_smiles_from_wurcs(wurcs_list, csdb_linear_cache, smiles_cache, glycoct_cache):
    sugar_smiles = {}
    updated_glycoct_cache = []
    updated_csdb_cache = []
    updated_smiles_cache = []
    for wurcs in wurcs_list:
        smiles = None
        glycoct = get_glycoct_from_wurcs(wurcs, glycoct_cache)
        updated_glycoct_cache.append({"WURCS": wurcs, "glycoct": glycoct})
        if not pd.isna(glycoct):
            csdb = get_csdb_from_glycoct(glycoct, csdb_linear_cache)
            updated_csdb_cache.append({"glycoct": glycoct, "csdb": csdb})
            if not pd.isna(csdb):
                smiles = get_smiles_from_csdb(csdb, smiles_cache)
                updated_smiles_cache.append({"csdb": csdb, "descriptor" : smiles})
        sugar_smiles[wurcs] = smiles
    updated_glycoct_cache_df = pd.concat([pd.DataFrame(updated_glycoct_cache, columns = ["WURCS", "glycoct"]), glycoct_cache]).drop_duplicates()
    updated_csdb_cache_df = pd.concat([pd.DataFrame(updated_csdb_cache, columns = ["glycoct", "csdb"]), csdb_linear_cache]).drop_duplicates()
    updated_smiles_cache_df = pd.concat([pd.DataFrame(updated_smiles_cache, columns = ["csdb", "descriptor"]), smiles_cache]).drop_duplicates()
    return sugar_smiles, updated_glycoct_cache_df, updated_csdb_cache_df, updated_smiles_cache_df

def get_chem_comp_descriptors(ccd_doc, comp_id_list):
    ligand_descriptors = {}
    for ligand in comp_id_list:
        lig_descriptor = None
        lig_block = ccd_doc.find_block(ligand)
        if lig_block is not None:
            lig_descriptors = pd.DataFrame(lig_block.find_mmcif_category("_pdbx_chem_comp_descriptor."), columns = ["comp_id", "type", "program", "program_version", "descriptor"])
            lig_descriptors["descriptor"] = lig_descriptors.descriptor.str.strip("\"|';").str.replace(r"\n$","", regex = True)
            lig_descriptors = lig_descriptors.loc[lig_descriptors.type == "SMILES"]
            PandasTools.AddMoleculeColumnToFrame(lig_descriptors, smilesCol='descriptor', molCol='pdb_ROMol')
            lig_descriptors = lig_descriptors.loc[lig_descriptors.pdb_ROMol.isna() == False]
            if len(lig_descriptors) == 0:
                lig_descriptor = None
            else: 
                #preference is to use openeye descriptors where available. if not, revert to the first smiles string able to be loaded into RDkit.
                preferred_row = lig_descriptors.loc[lig_descriptors.program.str.startswith("OpenEye")]
                if not preferred_row.empty:
                    lig_descriptor = preferred_row.iloc[0].descriptor
                else:
                    # Otherwise, select the first row with a SMILES string
                    lig_descriptor = lig_descriptors.iloc[0].descriptor
            ligand_descriptors[ligand] = lig_descriptor
        else:
            ligand_descriptors[ligand] = None
    return ligand_descriptors

def process_sifts_ec_map(sifts_ec_mapping_file, ec_records_file):
    sifts_chains_ec = sifts_ec_mapping_file.loc[sifts_ec_mapping_file.EC_NUMBER != "?"].copy()
    #the sifts mapping often has quotes in ec numbers because of extraction with gemmi it seems, strip these e.g. 2ex1
    sifts_chains_ec["EC_NUMBER"] = sifts_chains_ec["EC_NUMBER"].str.strip("'")
    sifts_chains_uniprot = sifts_ec_mapping_file.groupby(["PDB", "CHAIN"]).agg({"ACCESSION": set}).reset_index().copy()
    sifts_chains_uniprot["ACCESSION"] = sifts_chains_uniprot["ACCESSION"].apply(lambda x: "|".join(x)) #join the list of uniprot accessions with a pipe for downstream neo4j integration
    sifts_chains_uniprot.rename(columns = {"ACCESSION" : "uniprot_accession"}, inplace = True)

    sifts_chains_ec = sifts_chains_ec[["PDB", "CHAIN", "EC_NUMBER"]].groupby(["PDB", "CHAIN"]).agg({"EC_NUMBER": set}).reset_index() #group these into a set of pdb chain associated ec's
    sifts_chains_ec["EC_NUMBER"] = sifts_chains_ec["EC_NUMBER"].apply(lambda x: ",".join(x)) #join the list of EC numbers into a single string for ec records function
    sifts_chains_ec["EC_NUMBER"] = sifts_chains_ec["EC_NUMBER"].str.strip(",") #strip any leading or trailing commas - eventually drop the "" for np.nan in df before doing join so this isnt problem
    sifts_chains_ec = get_updated_enzyme_records(sifts_chains_ec, ec_records_file, ec_col = "EC_NUMBER")
    sifts_chains_ec.rename(columns = {"EC_NUMBER": "protein_entity_ec"}, inplace = True)

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
    parser.add_argument('--scop2_domains_info_file', type=str, help='scop2 domains info file')
    parser.add_argument('--scop2_descriptions_file', type=str, help='scop2 descriptions file')
    parser.add_argument('--interpro_xml', type=str, help='interpro xml file')
    parser.add_argument('--cath_names', type=str, help='cath names file')
    parser.add_argument('--cddf', type=str, help='cath domain definition file')
    parser.add_argument('--glycoct_cache', type=str, help='glycoct cache file')
    parser.add_argument('--smiles_cache', type=str, help='smiles cache file')
    parser.add_argument('--csdb_linear_cache', type=str, help='csdb linear cache file')
    parser.add_argument('--enzyme_dat_file', type=str, help='enzyme.dat file')
    parser.add_argument('--enzyme_class_file', type=str, help='enzclass.txt file')
    parser.add_argument('--sifts_ec_mapping', type=str, help='sifts ec mapping file')

    args = parser.parse_args()  
    
    if args.glycoct_cache:
        glycoct_cache = pd.read_pickle(args.glycoct_cache)
        #sometimes a second nan value is added to the cache, this is a sanity check to remove these
        glycoct_cache[~(glycoct_cache['WURCS'].duplicated(keep=False) & glycoct_cache['glycoct'].isna())]
    else:
        glycoct_cache = pd.DataFrame(columns = ["WURCS", "glycoct"])
    if args.smiles_cache:
        smiles_cache = pd.read_pickle(args.smiles_cache)
        smiles_cache = smiles_cache.loc[~(smiles_cache.csdb.duplicated(keep = False) & smiles_cache.descriptor.isna())]
    else:
        smiles_cache = pd.DataFrame(columns = ["csdb", "descriptor"])
    if args.csdb_linear_cache:
        csdb_linear_cache = pd.read_pickle(args.csdb_linear_cache)
        csdb_linear_cache.loc[~(csdb_linear_cache['glycoct'].duplicated(keep=False) & csdb_linear_cache['csdb'].isna())]  
    else:
        csdb_linear_cache = pd.DataFrame(columns = ["glycoct", "csdb"])

    contacts = pd.read_csv(args.contacts_file, sep = "\t", na_values = ["NaN", "None"], keep_default_na = False)   #for occurences of sodium

    ccd_doc = cif.read(args.ccd_cif)
    
    ligand_ids = contacts.loc[contacts.type == "ligand"].hetCode.unique()
    ligand_descriptors = get_chem_comp_descriptors(ccd_doc, ligand_ids)
    contacts.loc[contacts.type == "ligand", "descriptor"] = contacts.loc[contacts.type == "ligand"]["hetCode"].apply(lambda x: ligand_descriptors[x])

    wurcs_list = contacts.loc[contacts.type == "sugar", "descriptor"].unique()
    sugar_smiles, updated_glycoct_cache, updated_csdb_cache, updated_smiles_cache = get_sugar_smiles_from_wurcs(wurcs_list, csdb_linear_cache, smiles_cache, glycoct_cache)

    updated_glycoct_cache.to_pickle(f"glycoct_cache.pkl")
    updated_csdb_cache.to_pickle(f"csdb_cache.pkl")
    updated_smiles_cache.to_pickle(f"smiles_cache.pkl")

    contacts.loc[contacts.type == "sugar", "descriptor"] = contacts.loc[contacts.type == "sugar"]["descriptor"].apply(lambda x: sugar_smiles[x])
    #after merging the contacts with descriptors we probably want to flag or remove the ones with failed descriptors
    print(f"{len(contacts.loc[contacts.descriptor.isna() == True, 'uniqueID'].unique())} bound entities have failed to get a descriptor, removing")
    contacts = contacts.loc[contacts.descriptor.isna() == False].copy()
    ##then we need to assign a ligand uniqueID to all unique ligands in the contacts file - check implementation of this
    
    #assign ec and uniprot information to contacts

    ec_records_df_grouped = process_ec_records(args.enzyme_dat_file, args.enzyme_class_file)
    ec_records_df = ec_records_df_grouped.explode("ID")

    sifts_chains = pd.read_csv(f"{args.sifts_ec_mapping}", sep = "\t", comment="#")
    sifts_chains_ec, sifts_chains_uniprot = process_sifts_ec_map(sifts_chains, ec_records_df)
    contacts_ec = contacts.merge(sifts_chains_ec, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left", indicator = True)
    contacts_ec_unmatched = contacts_ec.loc[contacts_ec._merge != "both"].copy().drop(columns = ["_merge", "PDB", "CHAIN"])
    contacts_ec = contacts_ec.loc[contacts_ec._merge == "both"].copy().drop(columns = ["_merge", "PDB", "CHAIN"])
    contacts_ec_unmatched.to_csv(f"contacts_ec_unmatched.csv", index = False)

    #now add uniprot info to the contacts file
    contacts_ec_uniprot = contacts_ec.merge(sifts_chains_uniprot, left_on = ["pdb_id", "auth_chain_id"], right_on = ["PDB", "CHAIN"], how = "left", indicator = True)
    #what are we validating here to need the indicator col
    contacts_ec_uniprot.drop(columns = ["PDB", "CHAIN", "_merge"], inplace = True)
    
    #get the unique ligands from the dataset, and assign uniqueIDs - extract unique ligands to file for scoring.
    bound_entities_to_score = contacts_ec_uniprot[["descriptor", "description", "hetCode", "type", "ec_list"]].groupby(["hetCode", "description", "descriptor"]).agg({"ec_list": set}).reset_index().reset_index()
    #add ligand entity id to the contacts file
    contacts_ec_uniprot = contacts_ec_uniprot.merge(bound_entities_to_score[["descriptor", "description", "hetCode", "index"]], on = ["descriptor", "description", "hetCode"], how = "left", indicator = True)
    assert(len(contacts_ec_uniprot.loc[contacts_ec_uniprot._merge != "both"]) == 0)
    contacts_ec_uniprot.drop(columns = "_merge", inplace = True)
    contacts_ec_uniprot.rename(columns = {"index": "ligand_uniqueID"}, inplace = True)

    bound_entities_to_score.rename(columns = {"index" : "ligand_entity_id", "hetCode": "bl_name"}, inplace = True)
    bound_entities_to_score.to_pickle(f"bound_entities_to_score.pkl")

    # now that we have contacts_ec_uniprot, make a display EC list that aggregates all low level EC annotations into a range.
    

    #split the domain information into separate database specific dataframes for annotation
    cath_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "CATH"].copy()
    scop_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "SCOP"].copy()
    scop2b_fa_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "SCOP2B_Family"].copy()
    scop2b_sf_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "SCOP2B_SuperFamily"].copy()
    pfam_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "Pfam"].copy()
    gene3dsa_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "G3DSA"].copy()
    superfamily_contacts = contacts_ec_uniprot.loc[contacts_ec_uniprot.xref_db == "SuperFamily"].copy()

    core_cols = contacts_ec_uniprot.columns.tolist()
    cath_cols = ['cath_domain', 'cath_name', 'cath_code', 'cath_class', 'cath_class_name', 'cath_architecture', 'cath_architecture_name', 'cath_topology', 'cath_topology_name',
        'cath_homologous_superfamily', 'cath_homologous_superfamily_name']
    cath_mmcif_cols = ['cath_db_version', 'cath_db_verdate', 'cath_source', 'cath_domain_length', 'cath_domain_seq_header', 'cath_domain_seqs', 'cath_num_segments', 'cath_segments_dict']
    scop_cols = ['scop_id', 'sccs', 'domain_sunid', 'ancestor_sunid', 'cl_id', 'cf_id', 'sf_id', 'fa_id', 'dm_id', 'sp_id',
        'px_id', 'cl_description', 'cf_description', 'sf_description','fa_description', 'dm_description', 'sp_description', 'px_description']
    scop2b_fa_cols = ["SCOPCLA"]
    scop2b_sf_cols = ["SCOPCLA"]
    pfam_cols = ['clan', 'clan_acc', 'clan_comment', 'clan_description', 'pfam', 'pfam_accession', 'pfam_description', 'pfam_name']
    #interpro_cols = ['interpro_id','interpro_short_name', 'dbxref']
    #gene3dsa_cols
    #superfamily_cols
    #need to mkae the interpro cols be a part of the gene3dsa and superfamily cols - integrate on the derived from field

    #process domain database information
    if len(cath_contacts) > 0:
        cath_contacts_mmcif = cath_contacts.loc[cath_contacts.domain_type == "mmcif"].copy()
        if len(cath_contacts_mmcif) > 0:
            cath_domain_list = cath_contacts_mmcif.xref_db_acc.unique()
            cath_parsed_data = parse_cddf(args.cddf, cath_domain_list)
            cath_domains_info = build_cath_dataframe(cath_parsed_data)
            cath_contacts_mmcif = cath_contacts_mmcif.merge(cath_domains_info, how = "left", left_on = "xref_db_acc", right_on = "cath_domain", indicator = True)
            assert(len(cath_contacts_mmcif.loc[cath_contacts_mmcif._merge != "both"]) == 0)
            cath_contacts_mmcif.drop(columns = "_merge", inplace = True)
        else:
            cath_contacts_mmcif = pd.DataFrame(columns = core_cols + cath_cols + cath_mmcif_cols)
        cath_contacts_xml = cath_contacts.loc[cath_contacts.domain_type == "xml"].copy()
        if len(cath_contacts_xml) > 0:
            cath_domain_list = cath_contacts_xml.xref_db_acc.unique()
            cath_names = pd.read_csv(args.cath_names, sep = "    ", header = None, comment = "#", names = ["cath_code", "representative_domain", "name"])
            cath_names["name"] = cath_names["name"].str.replace("^:", "", regex = True)
            cath_domains_info = build_g3dsa_dataframe(cath_names,cath_domain_list)
            cath_contacts_xml = cath_contacts_xml.merge(cath_domains_info, how = "left", left_on = "xref_db_acc", right_on = "cath_domain", indicator = True)
            assert(len(cath_contacts.loc[cath_contacts_xml._merge != "both"]) == 0)
            cath_contacts_xml.drop(columns = "_merge", inplace = True)
            cath_contacts_xml[cath_mmcif_cols] = np.nan
        else:
            cath_contacts_xml = pd.DataFrame(columns = core_cols + cath_cols + cath_mmcif_cols)
        cath_contacts = pd.concat([cath_contacts_mmcif, cath_contacts_xml])
        cath_contacts.to_csv(f"cath_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        cath_contacts = pd.DataFrame(columns = core_cols + cath_cols + cath_mmcif_cols)
        cath_contacts.to_csv(f"cath_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")

    if len(scop_contacts) > 0:
        scop_domains_info = get_scop_domains_info(args.scop_domains_info_file, args.scop_descriptions_file)
        scop_contacts["xref_db_acc"] = scop_contacts["xref_db_acc"].astype(int)
        scop_contacts = scop_contacts.merge(scop_domains_info, how = "left", left_on = "xref_db_acc", right_on = "domain_sunid", indicator = True)
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

    if len(scop2b_fa_contacts) > 0:
        scop2b_fa_contacts.to_csv(f"scop2b_fa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        scop2b_fa_contacts = pd.DataFrame()
        scop2b_fa_contacts.to_csv(f"scop2b_fa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    scop2_domains_info = pd.read_csv(args.scop2_domains_info_file, sep = " ", comment = "#", header = None, names = ["FA-DOMID", "FA-PDBID","FA-PDBREG","FA-UNIID","FA-UNIREG","SF-DOMID","SF-PDBID","SF-PDBREG","SF-UNIID","SF-UNIREG","SCOPCLA"])
    scop2_domains_info[["SF-DOMID", "FA-DOMID"]] = scop2_domains_info[["SF-DOMID", "FA-DOMID"]].astype(str) #merge as string type with xref_db_acc in contacts 
    if len(scop2b_sf_contacts) > 0:
        scop2_sf_domains_info = scop2_domains_info[["SF-DOMID", "SCOPCLA"]].copy()
        scop2_sf_domains_info["SCOPCLA"] = scop2_sf_domains_info["SCOPCLA"].str.extract("(.*),FA=.*$") #remove the family level from sueprfamily level domains
        scop2_sf_domains_info = scop2_sf_domains_info.groupby("SF-DOMID").agg({"SCOPCLA": list}).reset_index()
        scop2_sf_domains_info["SCOPCLA"] = scop2_sf_domains_info["SCOPCLA"].str.join(";")

        scop2b_sf_contacts = scop2b_sf_contacts.merge(scop2_sf_domains_info, left_on = "xref_db_acc", right_on = "SF-DOMID", how = "left", indicator = True)
        #assert(len(scop2b_sf_contacts.loc[scop2b_sf_contacts._merge != "both"]) == 0) #we remoive the assertion here , because e.g. 8102391 from 8fvs is not present in the classification but is a domain. need to discuss this with SIFTS or SCOP team potentially
        scop2b_sf_contacts.drop(columns = ["_merge", "SF-DOMID"], inplace = True)
        scop2b_sf_contacts.to_csv(f"scop2b_sf_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        scop2b_sf_contacts = pd.DataFrame(columns = core_cols + scop2b_sf_cols)
        scop2b_sf_contacts.to_csv(f"scop2b_sf_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    if len(scop2b_fa_contacts) > 0:
        scop2_fa_domains_info = scop2_domains_info[["FA-DOMID", "SCOPCLA"]].copy()
        scop2_fa_domains_info = scop2_fa_domains_info.groupby("FA-DOMID").agg({"SCOPCLA": list}).reset_index()
        scop2_fa_domains_info["SCOPCLA"] = scop2_fa_domains_info["SCOPCLA"].str.join(";")

        scop2b_fa_contacts = scop2b_fa_contacts.merge(scop2_fa_domains_info, left_on = "xref_db_acc", right_on = "FA-DOMID", how = "left", indicator = True)
        #assert(len(scop2b_fa_contacts.loc[scop2b_fa_contacts._merge != "both"]) == 0) #as above 
        scop2b_fa_contacts.drop(columns = ["_merge", "FA-DOMID"], inplace = True)
        scop2b_fa_contacts.to_csv(f"scop2b_fa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        scop2b_fa_contacts = pd.DataFrame(columns = core_cols + scop2b_fa_cols)
        scop2b_fa_contacts.to_csv(f"scop2b_fa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    
    interpro_annotations = extract_interpro_domain_annotations(args.interpro_xml)
    
    if len(gene3dsa_contacts) > 0:
        #for gene3dsa we want to merge the
        gene3dsa_domain_list = gene3dsa_contacts.xref_db_acc.unique()
        cath_names = pd.read_csv(args.cath_names, sep = "    ", header = None, comment = "#", names = ["cath_code", "representative_domain", "name"])
        cath_names["name"] = cath_names["name"].str.replace("^:", "", regex = True)
        g3dsa_domains_info = build_g3dsa_dataframe(cath_names,gene3dsa_domain_list)
        gene3dsa_contacts = gene3dsa_contacts.merge(g3dsa_domains_info, how = "left", left_on = "xref_db_acc", right_on = "cath_domain", indicator = True)
        assert(len(gene3dsa_contacts.loc[gene3dsa_contacts._merge != "both"]) == 0)
        gene3dsa_contacts.drop(columns = "_merge", inplace = True) 
        gene3dsa_contacts = gene3dsa_contacts.merge(interpro_annotations, left_on = "derived_from", right_on = "interpro_accession", how = "left")
        gene3dsa_contacts.to_csv(f"gene3dsa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        gene3dsa_contacts = pd.DataFrame()
        gene3dsa_contacts.to_csv(f"gene3dsa_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    if len(superfamily_contacts) > 0:
        #needs to work on the scop1.75 data.
        superfamily_contacts = superfamily_contacts.merge(interpro_annotations, left_on = "derived_from", right_on = "interpro_accession", how = "left")
        superfamily_contacts.to_csv(f"superfamily_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")
    else:
        superfamily_contacts = pd.DataFrame()
        superfamily_contacts.to_csv(f"superfamily_pdb_residue_interactions.csv.gz", sep = "\t", index = False, compression = "gzip")

if __name__ == "__main__":
    main()