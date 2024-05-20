#!/usr/bin/env python

import pandas as pd
from rdkit import Chem
import numpy as np
import re
from Bio.ExPASy import Enzyme as EEnzyme
from pathlib import Path
import argparse
from utils import process_ec_records

import xml.etree.ElementTree as ET
import gzip

def sorted_set(x): 
    return [str(val) for val in sorted(set(map(int, x)))]

def get_all_chain_interactions(series):
    split_series = series.str.split("|")
    split_series = split_series.explode()
    split_series = split_series.astype(int)
    residues = [str(x) for x in sorted(split_series.unique())]
    residues = "|".join(residues)
    return residues

def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--enzyme_dat_file', metavar='enzyme_dat_file', type=str,
                        help = "path to enzyme.dat file")
    parser.add_argument('--enzyme_class_file', metavar='enzyme_class_file', type=str,
                        help = "path to enzclass.txt file")
    parser.add_argument('--cognate_ligands', metavar='cognate_ligands', type=str,
                        help = "path to biological ligands dataframe")
    parser.add_argument('--cath_domain_ownership', metavar='cath_domain_ownership', type=str,
                        help = "path to cath domain ownership file")
    parser.add_argument('--scop_domain_ownership', metavar='scop_domain_ownership', type=str,
                        help = "path to scop domain ownership file")
    parser.add_argument('--interpro_domain_ownership', metavar='interpro_domain_ownership', type=str,
                        help = "path to interpro domain ownership file")
    parser.add_argument('--pfam_domain_ownership', metavar='pfam_domain_ownership', type=str,
                        help = "path to pfam domain ownership file")
    parser.add_argument('--bound_entities', metavar='bound_entities', type=str, 
                        help = "path to bound entities to score df")
    parser.add_argument('--parity_calcs', metavar='parity_calcs', type=str,
                        help = "path to parity calcs dataframe")
    parser.add_argument('--parity_threshold', metavar='parity_threshold', type=float, default = 0.25,
                        help = "threshold for parity score")
    parser.add_argument('--rhea2ec', metavar='rhea2ec', type=str,
                        help = "path to rhea2ec file")
    parser.add_argument('--rhea_dir', metavar='rhea_dir', type=str,
                        help = "path to rhea directions file")
    parser.add_argument('--rhea_reaction_smiles', metavar='rhea_reaction_smiles', type=str,
                        help = "path to rhea reaction smiles file")

    args = parser.parse_args()

    ec_records_df_grouped = process_ec_records(args.enzyme_dat_file , args.enzyme_class_file)
    ec_id_nodes = ec_records_df_grouped[["TRANSFER", "DE"]].rename(columns = {"TRANSFER" : "ecID:ID(ec-id)", "DE" : "description"}).drop_duplicates()
    rhea2ec = pd.read_csv(f"{args.rhea2ec}", sep = "\t")
    rhea_dir = pd.read_csv(f"{args.rhea_dir}", sep = "\t")
    reaction_smiles = pd.read_csv(f"{args.rhea_reaction_smiles}", sep = "\t", header = None, names = ["RHEA_ID", "SMILES"])
    rheamerge = rhea2ec.merge(rhea_dir, left_on = "MASTER_ID" , right_on = "RHEA_ID_MASTER", how = "left")
    reactions_df_merged = reaction_smiles.merge(rheamerge[["RHEA_ID_LR", "ID"]], left_on = "RHEA_ID", right_on = "RHEA_ID_LR", how = "inner")
    reactions_df_merged["reactionSmiles"] = reactions_df_merged["RHEA_ID_LR"].astype("str") + ":" + reactions_df_merged["SMILES"]
    reactions_df_merged = reactions_df_merged.groupby("ID").agg({"reactionSmiles": list}).reset_index()
    reactions_df_merged["reactionSmiles"] = reactions_df_merged["reactionSmiles"].str.join("|")
    reactions_df_merged.rename(columns = {"reactionSmiles": "reactionSmiles:string[]"}, inplace = True)
    ec_id_nodes = ec_id_nodes.merge(reactions_df_merged[["ID", "reactionSmiles:string[]"]], left_on = "ecID:ID(ec-id)", right_on = "ID", how = "left", indicator = True)
    ec_nodes_class = ec_records_df_grouped[["class", "class_description"]].rename(columns = {"class": "ecID:ID(class-id)", "class_description": "description"}).drop_duplicates()
    ec_nodes_subclass = ec_records_df_grouped[["subclass", "subclass_description"]].rename(columns = {"subclass": "ecID:ID(subclass-id)", "subclass_description": "description"}).drop_duplicates()
    ec_nodes_subsubclass = ec_records_df_grouped[["subsubclass", "subsubclass_description"]].rename(columns = {"subsubclass": "ecID:ID(subsubclass-id)", "subsubclass_description": "description"}).drop_duplicates()
    ec_id_nodes.to_csv(f"ec_id_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_nodes_class.to_csv(f"ec_nodes_class.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_nodes_subclass.to_csv(f"ec_nodes_subclass.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_nodes_subsubclass.to_csv(f"ec_nodes_subsubclass.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_subsubclass_id_rel = ec_records_df_grouped[["TRANSFER", "subsubclass"]].drop_duplicates()
    ec_subsubclass_id_rel.rename(columns = {"TRANSFER" : ":START_ID(ec-id)", "subsubclass": ":END_ID(subsubclass-id)"}, inplace = True)
    ec_subclass_subsubclass_rel = ec_records_df_grouped[["subclass", "subsubclass"]].drop_duplicates()
    ec_subclass_subsubclass_rel.rename(columns = {"subsubclass" : ":START_ID(subsubclass-id)", "subclass": ":END_ID(subclass-id)"}, inplace = True)
    ec_class_subclass_rel = ec_records_df_grouped[["class", "subclass"]].drop_duplicates()
    ec_class_subclass_rel.rename(columns = {"subclass" : ":START_ID(subclass-id)", "class": ":END_ID(class-id)"}, inplace = True)
    ec_class_subclass_rel.to_csv(f"ec_class_subclass_rel.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_subclass_subsubclass_rel.to_csv(f"ec_subclass_subsubclass_rel.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_subsubclass_id_rel.to_csv(f"ec_subsubclass_id_rel.csv.gz", compression = "gzip", sep = "\t", index = False)

    
    cognate_ligands = pd.read_pickle(args.cognate_ligands)
    cognate_ligands_nodes = cognate_ligands[["canonical_smiles", "uniqueID", "ligand_db", "compound_name", "isCofactor", "compound_reaction"]].copy().drop_duplicates(subset = ["canonical_smiles", "uniqueID"])
    cognate_ligands_nodes.rename(columns = {"uniqueID": "uniqueID:ID(bio-id)", "compound_name": "name:string[]", "canonical_smiles": "canonicalSMILES", "ligand_db" : "ligandDB:string[]", "compound_reaction": "compoundReactionIDs:string[]"}, inplace = True)
    cognate_ligands_nodes.to_csv(f"cognate_ligand_nodes.csv.gz", compression = "gzip",sep = "\t", index = False)
    cognate_ligands_ec = cognate_ligands[["uniqueID", "entry"]].drop_duplicates()
    cognate_ligands_ec.rename(columns = {"uniqueID": ":START_ID(bio-id)", "entry": ":END_ID(ec-id)"}, inplace = True)
    cognate_ligands_ec.to_csv(f"cognate_ligands_ec.csv.gz", compression = "gzip", sep = "\t", index = False)


    cath_domains = pd.read_csv(f"{args.cath_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    scop_domains = pd.read_csv(f"{args.scop_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    pfam_domains = pd.read_csv(f"{args.pfam_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    interpro_domains = pd.read_csv(f"{args.interpro_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")

    pdb_nodes = pd.concat([cath_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords", "protein_entity_ec", "ec_list"]], scop_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords", "protein_entity_ec", "ec_list"]], pfam_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords", "protein_entity_ec", "ec_list"]], interpro_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords", "protein_entity_ec", "ec_list"]]]).drop_duplicates()
    pdb_nodes["pdb_keywords"] = pdb_nodes["pdb_keywords"].str.replace("\n", " ", regex = True)
    pdb_nodes["pdb_title"] = pdb_nodes["pdb_title"].str.replace("\n", " ", regex = True)
    pdb_nodes["pdb_descriptor"] = pdb_nodes["pdb_descriptor"].str.replace("\n", " ", regex = True)
    pdb_nodes["ec_list"] = pdb_nodes["ec_list"].str.replace(",", "|")
    pdb_nodes.loc[(pdb_nodes.protein_entity_ec.str.contains("-") == False) & (pdb_nodes.protein_entity_ec != pdb_nodes.ec_list), "updatedEC"] = "True"
    pdb_nodes.loc[(pdb_nodes.protein_entity_ec.str.contains("-")), "partialEC"] = "True"
    pdb_nodes["partialEC"] = pdb_nodes["partialEC"].fillna("False")
    pdb_nodes["updatedEC"] = pdb_nodes["updatedEC"].fillna("False")
    pdb_nodes.rename(columns = {"pdb_id": "pdbEntry:ID(pdb-id)", "pdb_title": "title", "pdb_descriptor": "description", "pdb_keywords": "keywords", "ec_list" : "ecList:string[]", "protein_entity_ec": "originalEC"}, inplace = True)
    pdb_nodes.to_csv(f"pdb_entry_nodes.csv.gz", sep='\t', compression = "gzip", index = False)

    pdb_ec_rels = pdb_nodes[["pdbEntry:ID(pdb-id)","ecList:string[]"]].copy()
    pdb_ec_rels["ecList:string[]"] = pdb_ec_rels["ecList:string[]"].str.split("|")
    pdb_ec_rels = pdb_ec_rels.explode("ecList:string[]")
    print(f"{len(pdb_ec_rels)} before filtering")
    pdb_ec_rels = pdb_ec_rels.loc[pdb_ec_rels["ecList:string[]"].isin(ec_id_nodes["ecID:ID(ec-id)"].unique())].copy()
    print(f"{len(pdb_ec_rels)} after filtering")
    pdb_ec_rels = pdb_ec_rels.loc[(pdb_ec_rels["ecList:string[]"] != "") & (pdb_ec_rels["ecList:string[]"].isna() == False)].reset_index()
    pdb_ec_rels.rename(columns = {"pdbEntry:ID(pdb-id)": ":START_ID(pdb-id)", "ecList:string[]": ":END_ID(ec-id)"}, inplace = True)
    pdb_ec_rels.to_csv(f"pdb_ec_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domains_nodes = scop_domains[["domain_accession", "assembly_chain_id_protein", "scop_id", "dm_description", "sccs", "domain_sunid"]].drop_duplicates()
    scop_domains_nodes["type"] = "SCOP"
    scop_domains_nodes["url"] = "https://scop.berkeley.edu/sunid=" + scop_domains_nodes["domain_sunid"].astype("str") + "&ver=1.75"
    scop_domains_nodes[":LABEL"] = scop_domains_nodes["type"] + "|domain"
    scop_domains_nodes.rename(columns = {"assembly_chain_id_protein": "assemblyChainID", "domain_accession": "domain:ID(scop-domain-id)", "scop_id": "scopAccession", "dm_description": "name", "sccs": "SCCS", "domain_sunid": "domainSUNID"}, inplace = True)
    scop_domains_nodes.to_csv(f"scop_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_domains_nodes = cath_domains[["domain_accession", "assembly_chain_id_protein", "cath_domain", "cath_name"]].drop_duplicates()
    cath_domains_nodes["type"] = "CATH"
    cath_domains_nodes["url"] = "https://www.cathdb.info/version/latest/superfamily/" + cath_domains_nodes["cath_domain"]
    cath_domains_nodes[":LABEL"] = cath_domains_nodes["type"] + "|domain"
    cath_domains_nodes.rename(columns = {"assembly_chain_id_protein": "assemblyChainID", "domain_accession": "domain:ID(cath-domain-id)", "cath_domain": "cathAccession", "cath_name": "name", "cath_homologous_superfamily": "homologousSuperfamily"}, inplace = True)
    cath_domains_nodes.to_csv(f"cath_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    interpro_domain_nodes = interpro_domains[["assembly_chain_id_protein", "domain_accession", 'interpro_accession','interpro_name']].drop_duplicates()
    interpro_domain_nodes["type"] = "InterProHomologousSuperfamily"
    interpro_domain_nodes["url"] = "https://www.ebi.ac.uk/interpro/entry/" + interpro_domain_nodes["interpro_accession"]
    interpro_domain_nodes[":LABEL"] = interpro_domain_nodes["type"] + "|domain"
    interpro_domain_nodes.rename(columns = {"assembly_chain_id_protein": "assemblyChainID", "domain_accession": "domain:ID(interpro-domain-id)" , "interpro_accession": "interproAccession", "interpro_name": "name"}, inplace = True)
    interpro_domain_nodes.to_csv(f"interpro_domain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    pfam_domains_nodes = pfam_domains[["assembly_chain_id_protein", "domain_accession", 'pfam_accession', 'pfam_name', 'pfam_description']].drop_duplicates()
    pfam_domains_nodes["type"] = "Pfam"
    pfam_domains_nodes["url"] = "https://www.ebi.ac.uk/interpro/entry/pfam/" + pfam_domains_nodes["pfam_accession"]
    pfam_domains_nodes[":LABEL"] = pfam_domains_nodes["type"] + "|domain"
    pfam_domains_nodes.rename(columns = {"assembly_chain_id_protein": "assemblyChainID", "domain_accession": "domain:ID(pfam-domain-id)", "pfam_accession": "pfamAccession", "pfam_description": "description", "pfam_name": "name"}, inplace = True)
    pfam_domains_nodes.to_csv(f"pfam_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    #domain_nodes = pd.concat([scop_domains_nodes, cath_domains_nodes, interpro_domain_nodes])
    #domain_nodes.to_csv(f"domain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_family_nodes = scop_domains[["sccs", "fa_id", "fa_description"]].drop_duplicates()
    scop_family_nodes.rename(columns = {"sccs": "SCCS", "fa_id": "scopFamily:ID(scop-family-id)", "fa_description": "description"}, inplace = True)
    scop_family_nodes.to_csv(f"scop_family_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_superfamily_nodes = scop_domains[["sccs", "sf_id", "sf_description"]]
    scop_superfamily_nodes["sccs"] = scop_superfamily_nodes["sccs"].str.extract(r"^(\w+\.\w+\.\w+)\.")
    scop_superfamily_nodes = scop_superfamily_nodes.drop_duplicates()
    scop_superfamily_nodes.rename(columns = {"sccs": "SCCS", "sf_id": "scopSuperfamily:ID(scop-superfam-id)", "sf_description": "description"}, inplace = True)
    scop_superfamily_nodes.to_csv(f"scop_superfamily_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_class_nodes = scop_domains[["sccs", "cl_id", "cl_description"]]
    scop_class_nodes["sccs"] = scop_class_nodes["sccs"].str.extract(r"(\w+)\.")
    scop_class_nodes = scop_class_nodes.drop_duplicates()
    scop_class_nodes.rename(columns = {"sccs": "SCCS", "cl_id": "scopClass:ID(scop-class-id)", "cl_description": "description"}, inplace = True)
    scop_class_nodes.to_csv(f"scop_class_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_fold_nodes = scop_domains[["sccs", "cf_id", "cf_description"]]
    scop_fold_nodes["sccs"] = scop_fold_nodes["sccs"].str.extract(r"^(\w+\.\w+)\.")
    scop_fold_nodes = scop_fold_nodes.drop_duplicates()
    scop_fold_nodes.rename(columns = {"sccs": "SCCS", "cf_id": "scopFold:ID(scop-fold-id)", "cf_description": "description"}, inplace = True)
    scop_fold_nodes.to_csv(f"scop_fold_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domain_family_rels = scop_domains[["domain_accession", "fa_id"]].drop_duplicates()
    scop_domain_family_rels.rename(columns = {"domain_accession": ":START_ID(scop-domain-id)", "fa_id": ":END_ID(scop-family-id)"}, inplace = True)
    scop_domain_family_rels.to_csv(f"scop_domain_family_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_family_superfamily_rels = scop_domains[["fa_id", "sf_id"]].drop_duplicates()
    scop_family_superfamily_rels.rename(columns = {"fa_id": ":START_ID(scop-family-id)", "sf_id": ":END_ID(scop-superfam-id)"}, inplace = True)
    scop_family_superfamily_rels.to_csv(f"scop_family_superfam_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_superfamily_fold_rels = scop_domains[["sf_id", "cf_id"]].drop_duplicates()
    scop_superfamily_fold_rels.rename(columns = {"sf_id": ":START_ID(scop-superfam-id)", "cf_id": ":END_ID(scop-fold-id)"}, inplace = True)
    scop_superfamily_fold_rels.to_csv(f"scop_superfam_fold_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_fold_class_rels = scop_domains[["cf_id", "cl_id"]].drop_duplicates()
    scop_fold_class_rels.rename(columns = {"cf_id": ":START_ID(scop-fold-id)", "cl_id": ":END_ID(scop-class-id)"}, inplace = True)
    scop_fold_class_rels.to_csv(f"scop_fold_class_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_class_nodes = cath_domains[["cath_class", "cath_class_name"]].drop_duplicates()
    cath_class_nodes.rename(columns = {"cath_class": "cathClass:ID(cath-class-ID)", "cath_class_name": "description"}, inplace = True)
    cath_class_nodes.to_csv(f"cath_class_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_architecture_nodes = cath_domains[["cath_architecture", "cath_architecture_name"]].drop_duplicates()
    cath_architecture_nodes.rename(columns = {"cath_architecture": "cathArchitecture:ID(cath-architecture-ID)", "cath_architecture_name": "description"}, inplace = True)
    cath_architecture_nodes.to_csv(f"cath_architecture_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_topology_nodes = cath_domains[["cath_topology", "cath_topology_name"]].drop_duplicates()
    cath_topology_nodes.rename(columns = {"cath_topology": "cathTopology:ID(cath-topology-ID)", "cath_topology_name": "description"}, inplace = True)
    cath_topology_nodes.to_csv(f"cath_topology_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_homologous_superfamily_nodes = cath_domains[["cath_homologous_superfamily", "cath_homologous_superfamily_name"]].drop_duplicates()
    cath_homologous_superfamily_nodes.rename(columns = {"cath_homologous_superfamily": "cathHomologousSuperfamily:ID(cath-homologous-superfamily-ID)", "cath_homologous_superfamily_name": "description"}, inplace = True)
    cath_homologous_superfamily_nodes.to_csv(f"cath_homologous_superfamily_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_class_architecture_rels = cath_domains[["cath_class", "cath_architecture"]].rename(columns = {"cath_class": ":END_ID(cath-class-ID)", "cath_architecture" : ":START_ID(cath-architecture-ID)"}).drop_duplicates()
    cath_architecture_topology_rels = cath_domains[["cath_architecture", "cath_topology"]].rename(columns = {"cath_architecture": ":END_ID(cath-architecture-ID)", "cath_topology" : ":START_ID(cath-topology-ID)"}).drop_duplicates()
    cath_topology_homology_rels = cath_domains[["cath_topology", "cath_homologous_superfamily"]].rename(columns = {"cath_topology": ":END_ID(cath-topology-ID)", "cath_homologous_superfamily" : ":START_ID(cath-homologous-superfamily-ID)"}).drop_duplicates()
    cath_homologous_superfamily_domain_rels = cath_domains[["cath_homologous_superfamily", "domain_accession"]].rename(columns = {"domain_accession": ":START_ID(cath-domain-id)", "cath_homologous_superfamily" : ":END_ID(cath-homologous-superfamily-ID)"}).drop_duplicates()

    cath_class_architecture_rels.to_csv(f"cath_class_architecture_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_architecture_topology_rels.to_csv(f"cath_architecture_topology_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_topology_homology_rels.to_csv(f"cath_topology_homology_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_homologous_superfamily_domain_rels.to_csv(f"cath_homologous_superfamily_domain_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    pfam_clans = pfam_domains.loc[(pfam_domains.clan_acc.isna() == False) & (pfam_domains.clan_acc != ""), ["clan_acc", "clan_description", "clan_comment"]].drop_duplicates()
    pfam_clans.rename(columns = {"clan_acc": "clanID:ID(pfam-clan-id)", "clan_description": "name", "clan_comment": "description"}, inplace = True)
    pfam_clans.to_csv(f"pfam_clans.csv.gz", compression = "gzip", sep = "\t", index = False)
    pfam_clan_rels = pfam_domains.loc[(pfam_domains.clan_acc.isna() == False) & (pfam_domains.clan_acc != ""), ["domain_accession", "clan_acc"]].rename(columns = {"domain_accession": ":START_ID(pfam-domain-id)", "clan_acc" : ":END_ID(pfam-clan-id)"}).drop_duplicates()
    pfam_clan_rels.to_csv(f"pfam_clan_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    bound_entities = pd.concat([cath_domains[["bound_ligand_struct_asym_id", "assembly_chain_id_ligand", 'bound_entity_pdb_residues', 'bound_molecule_display_id', 'hetCode', 'description', 'uniqueID', 'ligand_uniqueID', 'type' ,"ec_list"]],
                            scop_domains[["bound_ligand_struct_asym_id", "assembly_chain_id_ligand",'bound_entity_pdb_residues', 'bound_molecule_display_id', 'hetCode', 'description', 'uniqueID', 'ligand_uniqueID' , 'type' , "ec_list"]], 
                            pfam_domains[["bound_ligand_struct_asym_id", "assembly_chain_id_ligand",'bound_entity_pdb_residues', 'bound_molecule_display_id', 'hetCode', 'description', 'uniqueID', 'ligand_uniqueID', 'type' , "ec_list"]], 
                                interpro_domains[["bound_ligand_struct_asym_id", "assembly_chain_id_ligand", 'bound_entity_pdb_residues', 'bound_molecule_display_id', 'hetCode', 'description', 'uniqueID', 'ligand_uniqueID', 'type', "ec_list"]]]).drop_duplicates() #'represents', 

    #joining multiple ec lists in the instance where two chains with different ec annotations interact with a ligand.
    bound_entities = bound_entities.groupby([col for col in bound_entities.columns if col not in ["ec_list"]]).agg({"ec_list": set}).reset_index() 
    bound_entities["ec_list"] = bound_entities.ec_list.str.join("|")
    bound_entities["ec_list"] = bound_entities.ec_list.str.replace(",", "|") #replace comma separated lists that already existed with pipe delimiter

    bound_entities.rename(columns = {"uniqueID": "uniqueID:ID(be-id)", "ec_list": "ecList:string[]", "bound_entity_pdb_residues": "boundLigandResidues", "bound_ligand_struct_asym_id": "boundLigandStructChain", "assembly_chain_id_ligand": "boundLigandChain"}, inplace = True)
    bound_entities.loc[bound_entities.type == "sugar", "hetCode"] = "SUGAR"
    bound_entities["displayID"] = bound_entities["bound_molecule_display_id"] + ":" + bound_entities["hetCode"] + ":" + bound_entities["boundLigandResidues"].astype("str") + ":" + bound_entities["boundLigandChain"]

    bound_entities_descriptors = pd.read_pickle(args.bound_entities)
    bound_entities_descriptors["ec_list"] = bound_entities_descriptors.ec_list.str.join("|").str.replace(",", "|")

    bound_entities.rename(columns = {"ligand_uniqueID": "ligandUniqueID", "represents": "represents:string[]", "bound_molecule_display_id": "boundMoleculeID"}, inplace = True)
    
    bound_entities_descriptors.rename(columns = {"bl_name": "name", "ligand_entity_id": "uniqueID:ID(bd-id)"}, inplace = True)
    bound_entities_descriptors.drop(columns = "ec_list", inplace = True)
    bound_entities_descriptors.to_csv(f"bound_descriptors.csv.gz", compression = "gzip", sep = "\t", index = False)

    bound_entities_descriptors_rels = bound_entities[["uniqueID:ID(be-id)", "ligandUniqueID"]].rename(columns = {"uniqueID:ID(be-id)": ":START_ID(be-id)", "ligandUniqueID": ":END_ID(bd-id)"}).drop_duplicates()
    bound_entities_descriptors_rels.to_csv(f"be_bd_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    #bound_entities.drop(columns = ["ligandUniqueID"], inplace = True)
    bound_entities.to_csv(f"bound_entities.csv.gz", compression = "gzip", sep = "\t", index = False)

    bound_entities_pdb_rels = pd.concat([cath_domains[["uniqueID", "pdb_id"]], scop_domains[["uniqueID", "pdb_id"]], pfam_domains[["uniqueID", "pdb_id"]], interpro_domains[["uniqueID", "pdb_id"]]]).drop_duplicates()
    bound_entities_pdb_rels.rename(columns = {"uniqueID": ":START_ID(be-id)", "pdb_id": ":END_ID(pdb-id)"}, inplace = True)
    bound_entities_pdb_rels.to_csv(f"be_pdb_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    parity_calcs = pd.read_pickle(f"{args.parity_calcs}")
    parity_calcs["ec"] = parity_calcs["ec"].str.split(",")
    parity_calcs = parity_calcs.explode("ec")
    parity_calcs = parity_calcs.loc[parity_calcs.error.isna()]
    parity_calcs["parity_match_pdb"] = parity_calcs.parity_match.apply(lambda x: [str(x) for x in x.keys()] if isinstance(x, dict) else [])
    parity_calcs["parity_match_cognate"] = parity_calcs.parity_match.apply(lambda x: [str(x) for x in x.values()] if isinstance(x, dict) else [])
    parity_calcs["parity_match_pdb"] = parity_calcs["parity_match_pdb"].str.join("|")
    parity_calcs["parity_match_cognate"] = parity_calcs["parity_match_cognate"].str.join("|")

    parity_calcs_filtered = parity_calcs.loc[parity_calcs.score >= args.parity_threshold]

    parity_rels = bound_entities[["uniqueID:ID(be-id)", "ligandUniqueID", "ecList:string[]"]].copy()
    parity_rels["ecList:string[]"] = parity_rels["ecList:string[]"].str.split("|")
    parity_rels = parity_rels.explode("ecList:string[]")
    parity_rels = parity_rels.merge(parity_calcs_filtered, left_on = ["ligandUniqueID", "ecList:string[]"], right_on = ["pdb_ligand", "ec"], how = "inner")
    parity_rels = parity_rels[["uniqueID:ID(be-id)", "cognate_ligand", "pdbl_subparity", "score", "parity_smarts", "ec", "parity_match_pdb", "parity_match_cognate"]].rename(columns = {"uniqueID:ID(be-id)": ":START_ID(be-id)", "cognate_ligand": ":END_ID(bio-id)", "parity_match_pdb": "parityMatchPDB:int[]", "parity_match_cognate": "parityMatchCognate:int[]", "parity_smarts": "paritySMARTS"})
    parity_rels = parity_rels.groupby([col for col in parity_rels.columns if col != "ec"]).agg({"ec": list}).reset_index()
    parity_rels["ec"] = parity_rels["ec"].str.join("|")
    parity_rels["max_parity"] = parity_rels.groupby(":START_ID(be-id)")["score"].transform("max")
    parity_rels.loc[parity_rels.max_parity == parity_rels["score"], "bestCognate"] = "Y"
    parity_rels.bestCognate.fillna("N", inplace = True)
    parity_rels.drop(columns = ["max_parity"], inplace = True)
    parity_rels.rename(columns = {"score": "parityScore:float", "pdbl_subparity": "subParityScore:float", "ec": "ecList:string[]"}, inplace = True)
    parity_rels.to_csv(f"bound_entity_parity_score_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    if len(cath_domains) > 0:
        cath_domain_ligand_interactions = cath_domains[["assembly_chain_id_protein", "domain_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID", "bound_ligand_residue_interactions","domain_residue_interactions"]].drop_duplicates()
        cath_interface_mapping = cath_domain_ligand_interactions[["uniqueID", "assembly_chain_id_protein", "domain_residue_interactions"]].copy()
        cath_interface_mapping["allProteinInterface"] = cath_interface_mapping["domain_residue_interactions"].astype("str").str.split("|")
        cath_interface_mapping["allProteinInterface"] = cath_interface_mapping.apply(lambda x: "|".join([x.assembly_chain_id_protein + ":" + y for y in x.allProteinInterface]), axis = 1)
        cath_interface_mapping = cath_interface_mapping.groupby("uniqueID").agg({"allProteinInterface" : list}).reset_index()
        cath_interface_mapping["allProteinInterface"] = cath_interface_mapping.allProteinInterface.str.join("|")
        cath_domain_ligand_interactions = cath_domain_ligand_interactions.merge(cath_interface_mapping, how = "left", on = "uniqueID")
        cath_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "domain_accession": ":START_ID(cath-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_residue_interactions": "ligandInterface:string[]", "domain_residue_interactions": "proteinInterface:string[]", "allProteinInterface" : "allProteinInterface:string[]"}, inplace = True)
    else:
        cath_domain_ligand_interactions = pd.DataFrame(columns = [":END_ID(be-id)", ":START_ID(cath-domain-id)"])
    cath_domain_ligand_interactions.to_csv(f"cath_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)
    
    if len(scop_domains) > 0:
        scop_domain_ligand_interactions = scop_domains[["assembly_chain_id_protein", "domain_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID","bound_ligand_residue_interactions","domain_residue_interactions"]].drop_duplicates()
        scop_interface_mapping = scop_domain_ligand_interactions[["uniqueID", "assembly_chain_id_protein", "domain_residue_interactions"]].copy()
        scop_interface_mapping["allProteinInterface"] = scop_interface_mapping["domain_residue_interactions"].astype("str").str.split("|")
        scop_interface_mapping["allProteinInterface"] = scop_interface_mapping.apply(lambda x: "|".join([x.assembly_chain_id_protein + ":" + y for y in x.allProteinInterface]), axis = 1)
        scop_interface_mapping = scop_interface_mapping.groupby("uniqueID").agg({"allProteinInterface" : list}).reset_index()
        scop_interface_mapping["allProteinInterface"] = scop_interface_mapping.allProteinInterface.str.join("|")
        scop_domain_ligand_interactions = scop_domain_ligand_interactions.merge(scop_interface_mapping, how = "left", on = "uniqueID")
        scop_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "domain_accession": ":START_ID(scop-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_residue_interactions": "ligandInterface:string[]", "domain_residue_interactions": "proteinInterface:string[]", "allProteinInterface" : "allProteinInterface:string[]"}, inplace = True)
    else:
        scop_domain_ligand_interactions = pd.DataFrame(columns = [":END_ID(be-id)", ":START_ID(scop-domain-id)"])
    scop_domain_ligand_interactions.to_csv(f"scop_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)
    
    if len(pfam_domains) > 0:
        pfam_domain_ligand_interactions = pfam_domains[["assembly_chain_id_protein", "domain_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID", "bound_ligand_residue_interactions","domain_residue_interactions"]].drop_duplicates()
        pfam_interface_mapping = pfam_domain_ligand_interactions[["uniqueID", "assembly_chain_id_protein", "domain_residue_interactions"]].copy()
        pfam_interface_mapping["allProteinInterface"] = pfam_interface_mapping["domain_residue_interactions"].astype("str").str.split("|")
        pfam_interface_mapping["allProteinInterface"] = pfam_interface_mapping.apply(lambda x: "|".join([x.assembly_chain_id_protein + ":" + y for y in x.allProteinInterface]), axis = 1)
        pfam_interface_mapping = pfam_interface_mapping.groupby("uniqueID").agg({"allProteinInterface" : list}).reset_index()
        pfam_interface_mapping["allProteinInterface"] = pfam_interface_mapping.allProteinInterface.str.join("|")
        pfam_domain_ligand_interactions = pfam_domain_ligand_interactions.merge(pfam_interface_mapping, how = "left", on = "uniqueID")
        pfam_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "domain_accession": ":START_ID(pfam-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_residue_interactions": "ligandInterface:string[]", "domain_residue_interactions": "proteinInterface:string[]", "allProteinInterface" : "allProteinInterface:string[]"}, inplace = True)
    else:
        pfam_domain_ligand_interactions = pd.DataFrame(columns = [":END_ID(be-id)", ":START_ID(pfam-domain-id)"])
    pfam_domain_ligand_interactions.to_csv(f"pfam_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)
    
    if len(interpro_domains) > 0:
        interpro_domain_ligand_interactions = interpro_domains[["assembly_chain_id_protein", "domain_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID", "bound_ligand_residue_interactions","domain_residue_interactions"]].drop_duplicates()
        interpro_interface_mapping = interpro_domain_ligand_interactions[["uniqueID", "assembly_chain_id_protein", "domain_residue_interactions"]].copy()
        interpro_interface_mapping["allProteinInterface"] = interpro_interface_mapping["domain_residue_interactions"].astype("str").str.split("|")
        interpro_interface_mapping["allProteinInterface"] = interpro_interface_mapping.apply(lambda x: "|".join([x.assembly_chain_id_protein + ":" + y for y in x.allProteinInterface]), axis = 1)
        interpro_interface_mapping = interpro_interface_mapping.groupby("uniqueID").agg({"allProteinInterface" : list}).reset_index()
        interpro_interface_mapping["allProteinInterface"] = interpro_interface_mapping.allProteinInterface.str.join("|")
        interpro_domain_ligand_interactions = interpro_domain_ligand_interactions.merge(interpro_interface_mapping, how = "left", on = "uniqueID")
        interpro_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "domain_accession": ":START_ID(interpro-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_residue_interactions": "ligandInterface:string[]", "domain_residue_interactions": "proteinInterface:string[]", "allProteinInterface" : "allProteinInterface:string[]"}, inplace = True)        
    else:
        interpro_domain_ligand_interactions = pd.DataFrame(columns = [":END_ID(be-id)",":START_ID(interpro-domain-id)"])
    interpro_domain_ligand_interactions.to_csv(f"interpro_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)
    
    cath_pdb_rels = cath_domains[["domain_accession", "pdb_id"]].drop_duplicates().rename(columns = {"domain_accession": ":START_ID(cath-domain-id)", "pdb_id":":END_ID(pdb-id)"})
    scop_pdb_rels = scop_domains[["domain_accession", "pdb_id"]].drop_duplicates().rename(columns = {"domain_accession": ":START_ID(scop-domain-id)", "pdb_id":":END_ID(pdb-id)"})
    pfam_pdb_rels = pfam_domains[["domain_accession", "pdb_id"]].drop_duplicates().rename(columns = {"domain_accession": ":START_ID(pfam-domain-id)", "pdb_id":":END_ID(pdb-id)"})
    interpro_pdb_rels = interpro_domains[["domain_accession", "pdb_id"]].drop_duplicates().rename(columns = {"domain_accession": ":START_ID(interpro-domain-id)", "pdb_id":":END_ID(pdb-id)"})

    cath_pdb_rels.to_csv(f"cath_pdb_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    scop_pdb_rels.to_csv(f"scop_pdb_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    pfam_pdb_rels.to_csv(f"pfam_pdb_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    interpro_pdb_rels.to_csv(f"interpro_pdb_rels.csv.gz", compression = "gzip", sep = "\t", index = False)


    procoggraph_node = pd.DataFrame({"procoggraph:ID(procoggraph-id)": ["procoggraph"],
                                    "name": ["ProCogGraph"],
                                    "description": ["procoggraph"],
                                    "date_created": ["2024"],
                                    "date_updated": ["2024"],
                                    "database_version": ["0.1"],
                                    "biological_ligands_version": ["0.1"],
                                    "pdbe_graph_version": ["0.1"],
                                    "pdbe_graph_scripts_version": ["0.1"],
                                    "pdbe_graph_data_version": ["0.1"],
                                    "input_params": ["-"],})

    procoggraph_node.to_csv(f"procoggraph_node.csv.gz", compression = "gzip", sep = "\t", index = False)

if __name__ == "__main__":
    main()