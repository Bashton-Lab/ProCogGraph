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
    return [str(val) for val in sorted(set(x))]

def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--enzyme_dat_file', metavar='enzyme_dat_file', type=str,
                        help = "path to enzyme.dat file")
    parser.add_argument('--enzyme_class_file', metavar='enzyme_class_file', type=str,
                        help = "path to enzclass.txt file")
    parser.add_argument('--outdir', metavar='outdir', type=str,
                        help = "path to output directory")
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
    parser.add_argument('--bound_ligand_descriptors', metavar='bound_ligand_descriptors', type=str, 
                        help = "path to bound ligand descriptors dataframe")
    parser.add_argument('--bound_molecules_sugars_smiles', metavar='bound_molecules_sugars_smiles', type=str,
                        help = "path to bound molecules sugars ec dataframe")
    parser.add_argument('--parity_calcs', metavar='parity_calcs', type=str,
                        help = "path to parity calcs dataframe")

    args = parser.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    ec_records_df_grouped = process_ec_records(args.enzyme_dat_file , args.enzyme_class_file)

    ec_id_nodes = ec_records_df_grouped[["TRANSFER", "DE"]].rename(columns = {"TRANSFER" : "ecID:ID(ec-id)", "DE" : "description"}).drop_duplicates()
    ec_nodes_class = ec_records_df_grouped[["class", "class_description"]].rename(columns = {"class": "ecID:ID(class-id)", "class_description": "description"}).drop_duplicates()
    ec_nodes_subclass = ec_records_df_grouped[["subclass", "subclass_description"]].rename(columns = {"subclass": "ecID:ID(subclass-id)", "subclass_description": "description"}).drop_duplicates()
    ec_nodes_subsubclass = ec_records_df_grouped[["subsubclass", "subsubclass_description"]].rename(columns = {"subsubclass": "ecID:ID(subsubclass-id)", "subsubclass_description": "description"}).drop_duplicates()

    ec_id_nodes.to_csv(f"{args.outdir}/ec_id_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_nodes_class.to_csv(f"{args.outdir}/ec_nodes_class.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_nodes_subclass.to_csv(f"{args.outdir}/ec_nodes_subclass.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_nodes_subsubclass.to_csv(f"{args.outdir}/ec_nodes_subsubclass.csv.gz", compression = "gzip", sep = "\t", index = False)

    ec_subsubclass_id_rel = ec_records_df_grouped[["TRANSFER", "subsubclass"]].drop_duplicates()
    ec_subsubclass_id_rel.rename(columns = {"TRANSFER" : ":START_ID(ec-id)", "subsubclass": ":END_ID(subsubclass-id)"}, inplace = True)
    ec_subclass_subsubclass_rel = ec_records_df_grouped[["subclass", "subsubclass"]].drop_duplicates()
    ec_subclass_subsubclass_rel.rename(columns = {"subsubclass" : ":START_ID(subsubclass-id)", "subclass": ":END_ID(subclass-id)"}, inplace = True)
    ec_class_subclass_rel = ec_records_df_grouped[["class", "subclass"]].drop_duplicates()
    ec_class_subclass_rel.rename(columns = {"subclass" : ":START_ID(subclass-id)", "class": ":END_ID(class-id)"}, inplace = True)

    ec_class_subclass_rel.to_csv(f"{args.outdir}/ec_class_subclass_rel.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_subclass_subsubclass_rel.to_csv(f"{args.outdir}/ec_subclass_subsubclass_rel.csv.gz", compression = "gzip", sep = "\t", index = False)
    ec_subsubclass_id_rel.to_csv(f"{args.outdir}/ec_subsubclass_id_rel.csv.gz", compression = "gzip", sep = "\t", index = False)

    #It is questionable whether 204 properties is a good idea for a node. how can this be refactored out into nodes in the graph?
    cognate_ligands = pd.read_pickle(args.cognate_ligands)
    cognate_ligands_nodes = cognate_ligands[["canonical_smiles", "uniqueID", "ligand_db", "compound_name"]].copy().drop_duplicates(subset = ["canonical_smiles", "uniqueID"])
    cognate_ligands_nodes.rename(columns = {"uniqueID": "uniqueID:ID(bio-id)", "compound_name": "name", "canonical_smiles": "canonicalSMILES", "ligand_db" : "ligandDB:string[]"}, inplace = True)
    cognate_ligands_nodes.to_csv(f"{args.outdir}/cognate_ligand_nodes.csv.gz", compression = "gzip",sep = "\t", index = False)

    cognate_ligands_ec = cognate_ligands[["uniqueID", "entry"]].drop_duplicates()
    cognate_ligands_ec.rename(columns = {"uniqueID": ":START_ID(bio-id)", "entry": ":END_ID(ec-id)"}, inplace = True)
    cognate_ligands_ec.to_csv(f"{args.outdir}/cognate_ligands_ec.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_domains = pd.read_csv(f"{args.cath_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False)
    scop_domains = pd.read_csv(f"{args.scop_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False)
    pfam_domains = pd.read_csv(f"{args.pfam_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False)
    interpro_domains = pd.read_csv(f"{args.interpro_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False)

    pdb_nodes = pd.concat([cath_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords"]], scop_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords"]], pfam_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords"]], interpro_domains[["pdb_id", "pdb_title", "pdb_descriptor", "pdb_keywords"]]]).drop_duplicates()
    pdb_nodes["pdb_keywords"] = pdb_nodes["pdb_keywords"].str.replace("\n", " ", regex = True)
    pdb_nodes["pdb_title"] = pdb_nodes["pdb_title"].str.replace("\n", " ", regex = True)
    pdb_nodes["pdb_descriptor"] = pdb_nodes["pdb_descriptor"].str.replace("\n", " ", regex = True)
    pdb_nodes.rename(columns = {"pdb_id": "pdbEntry:ID(pdb-id)", "pdb_title": "title", "pdb_descriptor": "description", "pdb_keywords": "keywords"}, inplace = True)
    pdb_nodes.to_csv(f"{args.outdir}/pdb_entry_nodes.csv.gz", sep='\t', compression = "gzip", index = False)

    cath_protein_entities = cath_domains[["chainUniqueID", "chain_id",  "protein_entity_ec", "ec_list"]]
    scop_protein_entities = scop_domains[["chainUniqueID", "chain_id","protein_entity_ec", "ec_list"]]
    pfam_protein_entities = pfam_domains[["chainUniqueID", "chain_id", "protein_entity_ec", "ec_list"]]
    interpro_protein_entities = interpro_domains[["chainUniqueID", "chain_id", "protein_entity_ec", "ec_list"]]

    protein_entities = pd.concat([cath_protein_entities, scop_protein_entities, pfam_protein_entities, interpro_protein_entities])
    protein_entities = protein_entities.drop_duplicates(subset = ["chainUniqueID", "protein_entity_ec"])

    protein_entities["ec_list"] = protein_entities["ec_list"].str.replace(",", "|")
    protein_entities.loc[(protein_entities.protein_entity_ec.str.contains("-") == False) & (protein_entities.protein_entity_ec != protein_entities.ec_list), "updatedEC"] = "True"
    protein_entities.loc[(protein_entities.protein_entity_ec.str.contains("-")), "partialEC"] = "True"
    protein_entities["partialEC"] = protein_entities["partialEC"].fillna("False")
    protein_entities["updatedEC"] = protein_entities["updatedEC"].fillna("False")
    protein_entities.rename(columns = {"chainUniqueID": "pdbProteinChain:ID(pdbp-id)", "ec_list" : "ecList:string[]", "protein_entity_ec": "originalEC", "chain_id" : "chainID"}, inplace = True)
    protein_entities.to_csv(f"{args.outdir}/pdb_protein_chain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    protein_ec_rels = protein_entities[["pdbProteinChain:ID(pdbp-id)","ecList:string[]"]].copy()
    protein_ec_rels["ecList:string[]"] = protein_ec_rels["ecList:string[]"].str.split("|")
    protein_ec_rels = protein_ec_rels.explode("ecList:string[]")
    protein_ec_rels = protein_ec_rels.loc[(protein_ec_rels["ecList:string[]"] != "") & (protein_ec_rels["ecList:string[]"].isna() == False)].reset_index()
    protein_ec_rels.rename(columns = {"pdbProteinChain:ID(pdbp-id)": ":START_ID(pdbp-id)", "ecList:string[]": ":END_ID(ec-id)"}, inplace = True)
    protein_ec_rels.to_csv(f"{args.outdir}/protein_ec_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domains_nodes = scop_domains[["scop_id", "dm_description", "scop_sccs"]].drop_duplicates()
    scop_domains_nodes["type"] = "SCOP"
    scop_domains_nodes.rename(columns = {"scop_id": "domain:ID(scop-domain-id)", "dm_description": "name", "scop_sccs": "SCCS"}, inplace = True)
    scop_domains_nodes.to_csv(f"{args.outdir}/scop_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_domains_nodes = cath_domains[["cath_domain", "cath_name", "cath_homologous_superfamily"]].drop_duplicates()
    cath_domains_nodes["type"] = "CATH"
    cath_domains_nodes.rename(columns = {"cath_domain": "domain:ID(cath-domain-id)", "cath_name": "name", "cath_homologous_superfamily": "homologousSuperfamily"}, inplace = True)
    cath_domains_nodes.to_csv(f"{args.outdir}/cath_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    interpro_domain_nodes = interpro_domains[["interpro_accession", "interpro_name", "interpro_type", "dbxref"]].drop_duplicates()
    interpro_domain_nodes["interpro_type"] = "InterProHomologousSuperfamily"
    interpro_domain_nodes.rename(columns = {"interpro_accession": "domain:ID(interpro-domain-id)", "interpro_name": "name", "interpro_type": "type", "dbxref": "dbxref:string[]"}, inplace = True)
    interpro_domain_nodes.to_csv(f"{args.outdir}/interpro_domain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    pfam_domains_nodes = pfam_domains[['pfam_accession', 'pfam_name', 'pfam_description']].drop_duplicates()
    pfam_domains_nodes["type"] = "Pfam"
    pfam_domains_nodes.rename(columns = {"pfam_accession": "domain:ID(pfam-domain-id)", "pfam_description": "description", "pfam_name": "name"}, inplace = True)
    pfam_domains_nodes.to_csv(f"{args.outdir}/pfam_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    #domain_nodes = pd.concat([scop_domains_nodes, cath_domains_nodes, interpro_domain_nodes])
    #domain_nodes.to_csv(f"{args.outdir}/domain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_family_nodes = scop_domains[["fa_id", "fa_description"]].drop_duplicates()
    scop_family_nodes.rename(columns = {"fa_id": "scopFamily:ID(scop-family-id)", "fa_description": "description"}, inplace = True)
    scop_family_nodes.to_csv(f"{args.outdir}/scop_family_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_superfamily_nodes = scop_domains[["sf_id", "sf_description"]].drop_duplicates()
    scop_superfamily_nodes.rename(columns = {"sf_id": "scopSuperfamily:ID(scop-superfam-id)", "sf_description": "description"}, inplace = True)
    scop_superfamily_nodes.to_csv(f"{args.outdir}/scop_superfamily_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_class_nodes = scop_domains[["cl_id", "cl_description"]].drop_duplicates()
    scop_class_nodes.rename(columns = {"cl_id": "scopClass:ID(scop-class-id)", "cl_description": "description"}, inplace = True)
    scop_class_nodes.to_csv(f"{args.outdir}/scop_class_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_fold_nodes = scop_domains[["cf_id", "cf_description"]].drop_duplicates()
    scop_fold_nodes.rename(columns = {"cf_id": "scopFold:ID(scop-fold-id)", "cf_description": "description"}, inplace = True)
    scop_fold_nodes.to_csv(f"{args.outdir}/scop_fold_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domain_family_rels = scop_domains[["scop_id", "scop_sunid"]].drop_duplicates()
    scop_domain_family_rels.rename(columns = {"scop_id": ":START_ID(scop-domain-id)", "scop_sunid": ":END_ID(scop-family-id)"}, inplace = True)
    scop_domain_family_rels.to_csv(f"{args.outdir}/scop_domain_family_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_family_superfamily_rels = scop_domains[["scop_sunid", "sf_id"]].drop_duplicates()
    scop_family_superfamily_rels.rename(columns = {"scop_sunid": ":START_ID(scop-family-id)", "sf_id": ":END_ID(scop-superfam-id)"}, inplace = True)
    scop_family_superfamily_rels.to_csv(f"{args.outdir}/scop_family_superfam_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_superfamily_fold_rels = scop_domains[["sf_id", "cf_id"]].drop_duplicates()
    scop_superfamily_fold_rels.rename(columns = {"sf_id": ":START_ID(scop-superfam-id)", "cf_id": ":END_ID(scop-fold-id)"}, inplace = True)
    scop_superfamily_fold_rels.to_csv(f"{args.outdir}/scop_superfam_fold_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_fold_class_rels = scop_domains[["cf_id", "cl_id"]].drop_duplicates()
    scop_fold_class_rels.rename(columns = {"cf_id": ":START_ID(scop-fold-id)", "cl_id": ":END_ID(scop-class-id)"}, inplace = True)
    scop_fold_class_rels.to_csv(f"{args.outdir}/scop_fold_class_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_class_nodes = cath_domains.cath_class.unique()
    cath_architecture_nodes = cath_domains.cath_architecture.unique()
    cath_topology_nodes = cath_domains.cath_topology.unique()
    cath_homologous_superfamily_nodes = cath_domains.cath_homologous_superfamily.unique()

    np.savetxt(f"{args.outdir}/cath_class_nodes.csv.gz", cath_class_nodes, delimiter='\t',fmt='%s', header='cathClass:ID(cath-class-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_architecture_nodes.csv.gz", cath_architecture_nodes, delimiter='\t',fmt='%s', header='cathArchitecture:ID(cath-architecture-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_topology_nodes.csv.gz", cath_topology_nodes, delimiter='\t',fmt='%s', header='cathTopology:ID(cath-topology-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_homologous_superfamily_nodes.csv.gz", cath_homologous_superfamily_nodes, delimiter='\t',fmt='%s', header='cathHomologousSuperfamily:ID(cath-homologous-superfamily-ID)', comments='')

    cath_class_architecture_rels = cath_domains[["cath_class", "cath_architecture"]].rename(columns = {"cath_class": ":END_ID(cath-class-ID)", "cath_architecture" : ":START_ID(cath-architecture-ID)"}).drop_duplicates()
    cath_architecture_topology_rels = cath_domains[["cath_architecture", "cath_topology"]].rename(columns = {"cath_architecture": ":END_ID(cath-architecture-ID)", "cath_topology" : ":START_ID(cath-topology-ID)"}).drop_duplicates()
    cath_topology_homology_rels = cath_domains[["cath_topology", "cath_homologous_superfamily"]].rename(columns = {"cath_topology": ":END_ID(cath-topology-ID)", "cath_homologous_superfamily" : ":START_ID(cath-homologous-superfamily-ID)"}).drop_duplicates()
    cath_homologous_superfamily_domain_rels = cath_domains[["cath_homologous_superfamily", "cath_domain"]].rename(columns = {"cath_domain": ":START_ID(cath-domain-id)", "cath_homologous_superfamily" : ":END_ID(cath-homologous-superfamily-ID)"}).drop_duplicates()

    cath_class_architecture_rels.to_csv(f"{args.outdir}/cath_class_architecture_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_architecture_topology_rels.to_csv(f"{args.outdir}/cath_architecture_topology_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_topology_homology_rels.to_csv(f"{args.outdir}/cath_topology_homology_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_homologous_superfamily_domain_rels.to_csv(f"{args.outdir}/cath_homologous_superfamily_domain_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    pfam_clans = pfam_domains[["clan_acc", "clan_description", "clan_comment"]].drop_duplicates()
    pfam_clans.rename(columns = {"clan_acc": "clanID:ID(pfam-clan-id)", "clan_description": "name", "clan_comment": "description"}, inplace = True)
    pfam_clans.to_csv(f"{args.outdir}/pfam_clans.csv.gz", compression = "gzip", sep = "\t", index = False)
    pfam_clan_rels = pfam_domains.loc[pfam_domains.clan_acc.isna() == False, ["pfam_accession", "clan_acc"]].rename(columns = {"pfam_accession": ":START_ID(pfam-domain-id)", "clan_acc" : ":END_ID(pfam-clan-id)"}).drop_duplicates()
    pfam_clan_rels.to_csv(f"{args.outdir}/pfam_clan_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    bound_entities = pd.concat([cath_domains[['bm_ids', "bound_ligand_struct_asym_id", 'bound_ligand_auth_id', 'bound_molecule_display_id', 'name', 'description', 'uniqueID', 'ligand_uniqueID', 'bound_ligand_id',  'type' ,"ec_list"]],
                            scop_domains[['bm_ids', "bound_ligand_struct_asym_id",'bound_ligand_auth_id', 'bound_molecule_display_id', 'name', 'description', 'uniqueID', 'ligand_uniqueID', 'bound_ligand_id', 'type' , "ec_list"]], 
                            pfam_domains[['bm_ids', "bound_ligand_struct_asym_id",'bound_ligand_auth_id', 'bound_molecule_display_id', 'name', 'description', 'uniqueID', 'ligand_uniqueID', 'bound_ligand_id', 'type' , "ec_list"]], 
                                interpro_domains[['bm_ids', "bound_ligand_struct_asym_id", 'bound_ligand_auth_id', 'bound_molecule_display_id', 'name', 'description', 'uniqueID', 'ligand_uniqueID','bound_ligand_id',  'type', "ec_list"]]]).drop_duplicates()

    #sometimes a ligand is bound by two different protein chains with different ec. here we join the ec lists for these together.
    bound_entities["bound_ligand_auth_id"] = bound_entities["bound_ligand_auth_id"].astype("str").str.split("\|")
    bound_entities = bound_entities.explode("bound_ligand_auth_id")
    bound_entities["bound_ligand_auth_id"] = bound_entities["bound_ligand_auth_id"].astype("int")
    bound_entities = bound_entities.groupby([col for col in bound_entities.columns if col not in ["ec_list", "bound_ligand_id", "bound_ligand_auth_id", "bound_ligand_struct_asym_id"]]).agg({"ec_list": set, "bound_ligand_id": set, "bound_ligand_auth_id": sorted_set, "bound_ligand_struct_asym_id": "first"}).reset_index() #joining the multiple bound ligands making up a sugar, and multiple ec lists in the instance where two chains with different ec annotations interact with a ligand.
    bound_entities["bound_ligand_auth_id"] = bound_entities["bound_ligand_auth_id"].str.join("|")
    bound_entities["bound_ligand_id"] = bound_entities.bound_ligand_id.str.join("|")
    bound_entities["ec_list"] = bound_entities.ec_list.str.join("|")
    bound_entities["ec_list"] = bound_entities.ec_list.str.replace(",", "|") #replace comma separated lists that already existed with pipe delimiter


    bound_entities.rename(columns = {"uniqueID": "uniqueID:ID(be-id)", "name": "hetCode", "bound_ligand_id": "componentLigands:string[]", "ec_list": "ecList:string[]", "bound_ligand_auth_id": "boundLigandResidues", "bound_ligand_struct_asym_id": "boundLigandChain"}, inplace = True)
    bound_entities.loc[bound_entities.type == "sugar", "hetCode"] = "SUGAR"
    bound_entities["displayID"] = bound_entities["bound_molecule_display_id"] + ":" + bound_entities["hetCode"] + ":" + bound_entities["boundLigandResidues"] + ":" + bound_entities["boundLigandChain"]

    bound_ligand_descriptors = pd.read_pickle(args.bound_ligand_descriptors)
    bound_molecules_sugar_smiles = pd.read_pickle(args.bound_molecules_sugars_smiles)
    bound_sugar_descriptors = bound_molecules_sugar_smiles[["ligand_index", "description", "descriptor", "ec_list"]].copy()
    bound_sugar_descriptors = bound_sugar_descriptors.groupby(["ligand_index", "description", "descriptor"]).agg({"ec_list": set}).reset_index()
    bound_sugar_descriptors["bl_name"] = bound_sugar_descriptors["description"]
    bound_sugar_descriptors.rename(columns = {"ligand_index": "ligand_entity_id"}, inplace = True)
    bound_descriptors = pd.concat([bound_ligand_descriptors, bound_sugar_descriptors])
    bound_descriptors["ec_list"] = bound_descriptors.ec_list.str.join("|")

    bound_entities.rename(columns = {"ligand_uniqueID": "ligandUniqueID", "bm_ids": "boundMolecules:string[]", "bound_molecule_display_id": "boundMoleculeID"}, inplace = True)
    
    bound_descriptors.rename(columns = {"bl_name": "name", "ligand_entity_id": "uniqueID:ID(bd-id)"}, inplace = True)
    bound_descriptors.drop(columns = "ec_list", inplace = True)
    bound_descriptors.to_csv(f"{args.outdir}/bound_descriptors.csv.gz", compression = "gzip", sep = "\t", index = False)
    bound_entities_descriptors_rels = bound_entities[["uniqueID:ID(be-id)", "ligandUniqueID"]].rename(columns = {"uniqueID:ID(be-id)": ":START_ID(be-id)", "ligandUniqueID": ":END_ID(bd-id)"}).drop_duplicates()
    bound_entities_descriptors_rels.to_csv(f"{args.outdir}/be_bd_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    #bound_entities.drop(columns = ["ligandUniqueID"], inplace = True)
    bound_entities.to_csv(f"{args.outdir}/bound_entities.csv.gz", compression = "gzip", sep = "\t", index = False)

    bound_entities_pdb_rels = pd.concat([cath_domains[["uniqueID", "pdb_id"]], scop_domains[["uniqueID", "pdb_id"]], interpro_domains[["uniqueID", "pdb_id"]]]).drop_duplicates()
    bound_entities_pdb_rels.rename(columns = {"uniqueID": ":START_ID(be-id)", "pdb_id": ":END_ID(pdb-id)"}, inplace = True)
    bound_entities_pdb_rels.to_csv(f"{args.outdir}/be_pdb_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    parity_calcs = pd.read_pickle(f"{args.parity_calcs}")
    parity_calcs["ec"] = parity_calcs["ec"].str.split(",")
    parity_calcs = parity_calcs.explode("ec")
    parity_calcs_filtered = parity_calcs.loc[parity_calcs.score >= 0.25]

    parity_rels = bound_entities[["uniqueID:ID(be-id)", "ligandUniqueID", "ecList:string[]"]].copy()
    parity_rels["ecList:string[]"] = parity_rels["ecList:string[]"].str.split("|")
    parity_rels = parity_rels.explode("ecList:string[]")
    parity_rels = parity_rels.merge(parity_calcs_filtered, left_on = ["ligandUniqueID", "ecList:string[]"], right_on = ["pdb_ligand", "ec"], how = "inner")
    parity_rels = parity_rels.merge(cognate_ligands[["canonical_smiles", "uniqueID"]].drop_duplicates(), left_on = "compound", right_on = "canonical_smiles", how = "left", indicator = True)
    assert(len(parity_rels.loc[parity_rels._merge != "both"]) == 0)
    parity_rels.drop(columns = ["_merge"], inplace = True)
    parity_rels = parity_rels[["uniqueID:ID(be-id)", "uniqueID", "pdbl_subparity", "score", "ec"]].rename(columns = {"uniqueID:ID(be-id)": ":START_ID(be-id)", "uniqueID": ":END_ID(bio-id)"})
    parity_rels = parity_rels.groupby([col for col in parity_rels.columns if col != "ec"]).agg({"ec": list}).reset_index()
    parity_rels["ec"] = parity_rels["ec"].str.join("|")
    parity_rels.rename(columns = {"score": "parityScore:float", "pdbl_subparity": "subParityScore:float", "ec": "ecList:string[]"}, inplace = True)
    parity_rels.to_csv(f"{args.outdir}/bound_entity_parity_score_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_domain_ligand_interactions = cath_domains[["cath_domain", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID", "bound_ligand_auth_id","pdb_residue_auth_id"]].drop_duplicates()
    cath_domain_ligand_interactions["bound_ligand_auth_id"] = cath_domain_ligand_interactions["bound_ligand_auth_id"].astype("str").str.split("|")
    cath_domain_ligand_interactions["bound_ligand_auth_id"] = cath_domain_ligand_interactions["bound_ligand_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    cath_domain_ligand_interactions["pdb_residue_auth_id"] = cath_domain_ligand_interactions["pdb_residue_auth_id"].astype("str").str.split("|")
    cath_domain_ligand_interactions["pdb_residue_auth_id"] = cath_domain_ligand_interactions["pdb_residue_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    cath_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "cath_domain": ":START_ID(cath-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_auth_id": "ligandInterface", "pdb_residue_auth_id": "proteinInterface"}, inplace = True)
    cath_domain_ligand_interactions.to_csv(f"{args.outdir}/cath_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domain_ligand_interactions = scop_domains[["scop_id", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID","bound_ligand_auth_id","pdb_residue_auth_id"]].drop_duplicates()
    scop_domain_ligand_interactions["bound_ligand_auth_id"] = scop_domain_ligand_interactions["bound_ligand_auth_id"].astype("str").str.split("|")
    scop_domain_ligand_interactions["bound_ligand_auth_id"] = scop_domain_ligand_interactions["bound_ligand_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    scop_domain_ligand_interactions["pdb_residue_auth_id"] = scop_domain_ligand_interactions["pdb_residue_auth_id"].astype("str").str.split("|")
    scop_domain_ligand_interactions["pdb_residue_auth_id"] = scop_domain_ligand_interactions["pdb_residue_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    scop_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "scop_id": ":START_ID(scop-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_auth_id": "ligandInterface", "pdb_residue_auth_id": "proteinInterface"}, inplace = True)
    scop_domain_ligand_interactions.to_csv(f"{args.outdir}/scop_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)

    pfam_domain_ligand_interactions = pfam_domains[["pfam_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID", "bound_ligand_auth_id","pdb_residue_auth_id"]].drop_duplicates()
    pfam_domain_ligand_interactions["bound_ligand_auth_id"] = pfam_domain_ligand_interactions["bound_ligand_auth_id"].astype("str").str.split("|")
    pfam_domain_ligand_interactions["bound_ligand_auth_id"] = pfam_domain_ligand_interactions["bound_ligand_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    pfam_domain_ligand_interactions["pdb_residue_auth_id"] = pfam_domain_ligand_interactions["pdb_residue_auth_id"].astype("str").str.split("|")
    pfam_domain_ligand_interactions["pdb_residue_auth_id"] = pfam_domain_ligand_interactions["pdb_residue_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    pfam_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "pfam_accession": ":START_ID(pfam-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_auth_id": "ligandInterface", "pdb_residue_auth_id": "proteinInterface"}, inplace = True)
    pfam_domain_ligand_interactions.to_csv(f"{args.outdir}/pfam_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)

    interpro_domain_ligand_interactions = interpro_domains[["interpro_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_covalent_counts", "domain_ownership", "uniqueID", "bound_ligand_auth_id","pdb_residue_auth_id"]].drop_duplicates()
    interpro_domain_ligand_interactions["bound_ligand_auth_id"] = interpro_domain_ligand_interactions["bound_ligand_auth_id"].astype("str").str.split("|")
    interpro_domain_ligand_interactions["bound_ligand_auth_id"] = interpro_domain_ligand_interactions["bound_ligand_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    interpro_domain_ligand_interactions["pdb_residue_auth_id"] = interpro_domain_ligand_interactions["pdb_residue_auth_id"].astype("str").str.split("|")
    interpro_domain_ligand_interactions["pdb_residue_auth_id"] = interpro_domain_ligand_interactions["pdb_residue_auth_id"].apply(lambda x: sorted_set(x)).str.join("|")
    interpro_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "interpro_accession": ":START_ID(interpro-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_covalent_counts": "domainCovalentCounts", "domain_ownership" : "interactionMode", "bound_ligand_auth_id": "ligandInterface", "pdb_residue_auth_id": "proteinInterface"}, inplace = True)
    interpro_domain_ligand_interactions.to_csv(f"{args.outdir}/interpro_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_pdb_protein_rels = scop_domains[["pdb_id", "chainUniqueID"]].drop_duplicates()
    cath_pdb_protein_rels = cath_domains[["pdb_id", "chainUniqueID"]].drop_duplicates()
    pfam_pdb_protein_rels = pfam_domains[["pdb_id", "chainUniqueID"]].drop_duplicates()
    interpro_pdb_protein_rels = interpro_domains[["pdb_id", "chainUniqueID"]].drop_duplicates()

    pdb_protein_rels = pd.concat([scop_pdb_protein_rels, cath_pdb_protein_rels, pfam_pdb_protein_rels, interpro_pdb_protein_rels]).drop_duplicates()
    pdb_protein_rels.rename(columns = {"chainUniqueID": ":START_ID(pdbp-id)", "pdb_id": ":END_ID(pdb-id)"}, inplace = True)
    pdb_protein_rels.to_csv(f"{args.outdir}/pdb_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_protein_rels = cath_domains[["cath_domain", "chainUniqueID"]].drop_duplicates().rename(columns = {"cath_domain": ":START_ID(cath-domain-id)", "chainUniqueID":":END_ID(pdbp-id)"})
    scop_protein_rels = scop_domains[["scop_id", "chainUniqueID"]].drop_duplicates().rename(columns = {"scop_id": ":START_ID(scop-domain-id)", "chainUniqueID":":END_ID(pdbp-id)"})
    pfam_protein_rels = pfam_domains[["pfam_accession", "chainUniqueID"]].drop_duplicates().rename(columns = {"pfam_accession": ":START_ID(pfam-domain-id)", "chainUniqueID":":END_ID(pdbp-id)"})
    interpro_protein_rels = interpro_domains[["interpro_accession", "chainUniqueID"]].drop_duplicates().rename(columns = {"interpro_accession": ":START_ID(interpro-domain-id)", "chainUniqueID":":END_ID(pdbp-id)"})

    cath_protein_rels.to_csv(f"{args.outdir}/cath_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    scop_protein_rels.to_csv(f"{args.outdir}/scop_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    pfam_protein_rels.to_csv(f"{args.outdir}/pfam_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    interpro_protein_rels.to_csv(f"{args.outdir}/interpro_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

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

    procoggraph_node.to_csv(f"{args.outdir}/procoggraph_node.csv.gz", compression = "gzip", sep = "\t", index = False)

if __name__ == "__main__":
    main()