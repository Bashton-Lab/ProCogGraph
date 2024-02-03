#!/usr/bin/env python

import pandas as pd
from rdkit import Chem
import numpy as np
import re
from Bio.ExPASy import Enzyme as EEnzyme
from pathlib import Path
import argparse
from utils import process_ec_records

def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--enzyme_dat_file', metavar='enzyme_dat_file', type=str,
                        help = "path to enzyme.dat file")
    parser.add_argument('--enzyme_class_file', metavar='enzyme_class_file', type=str,
                        help = "path to enzclass.txt file")
    parser.add_argument('--outdir', metavar='outdir', type=str,
                        help = "path to output directory")
    parser.add_argument('--biological_ligands', metavar='biological_ligands', type=str,
                        help = "path to biological ligands dataframe")
    parser.add_argument('--cath_domain_ownership', metavar='cath_domain_ownership', type=str,
                        help = "path to cath domain ownership file")
    parser.add_argument('--scop_domain_ownership', metavar='scop_domain_ownership', type=str,
                        help = "path to scop domain ownership file")
    parser.add_argument('--interpro_domain_ownership', metavar='interpro_domain_ownership', type=str,
                        help = "path to interpro domain ownership file")
    parser.add_argument('--bound_ligand_descriptors', metavar='bound_ligand_descriptors', type=str, 
                        help = "path to bound ligand descriptors dataframe")
    parser.add_argument('--bound_molecules_sugars_smiles', metavar='bound_molecules_sugars_smiles', type=str,
                        help = "path to bound molecules sugars ec dataframe")
    parser.add_argument('--parity_calcs', metavar='parity_calcs', type=str,
                        help = "path to parity calcs dataframe")

    args = parser.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    #check if all files to be produced already exist in dir, exit if they do
    #if 

    ec_records_df_grouped = process_ec_records(args.enzyme_dat_file , args.enzyme_class_file)

    ec_id_nodes = ec_records_df_grouped[["TRANSFER", "DE"]].rename(columns = {"TRANSFER" : "ecID:ID(ec-id)", "DE" : "description"}).drop_duplicates()
    ec_nodes_class = ec_records_df_grouped[["class", "class_description"]].rename(columns = {"class": "ecClass:ID(class-id)", "class_description": "description"}).drop_duplicates()
    ec_nodes_subclass = ec_records_df_grouped[["subclass", "subclass_description"]].rename(columns = {"subclass": "ecSubClass:ID(subclass-id)", "subclass_description": "description"}).drop_duplicates()
    ec_nodes_subsubclass = ec_records_df_grouped[["subsubclass", "subsubclass_description"]].rename(columns = {"subsubclass": "ecSubSubClass:ID(subsubclass-id)", "subsubclass_description": "description"}).drop_duplicates()

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
    biological_ligands = pd.read_pickle(args.biological_ligands)

    biological_ligands_nodes = biological_ligands[[col for col in biological_ligands.columns if col not in ["entry", "_merge"]]].copy()
    biological_ligands_nodes.drop_duplicates(subset = ["canonical_smiles", "uniqueID"], inplace = True)
    biological_ligands_nodes = biological_ligands_nodes[["canonical_smiles", "uniqueID", "ligand_db", "compound_name"]]
    biological_ligands_nodes.rename(columns = {"uniqueID": "uniqueID:ID(bio-id)", "compound_name": "name", "canonical_smiles": "canonicalSMILES", "ligand_db" : "ligandDB:string[]"}, inplace = True)

    biological_ligands_nodes.to_csv(f"{args.outdir}/biological_ligand_nodes.csv.gz", compression = "gzip",sep = "\t", index = False)

    biological_ligands_ec = biological_ligands[["uniqueID", "entry"]].drop_duplicates()
    biological_ligands_ec.rename(columns = {"uniqueID": ":START_ID(bio-id)", "entry": ":END_ID(ec-id)"}, inplace = True)
    biological_ligands_ec.to_csv(f"{args.outdir}/biological_ligands_ec.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_domains = pd.read_csv(f"{args.cath_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False)
    scop_domains = pd.read_csv(f"{args.scop_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False)
    interpro_domains = pd.read_csv(f"{args.interpro_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False)

    pdb_nodes = pd.concat([cath_domains.pdb_id, scop_domains.pdb_id, interpro_domains.pdb_id]).unique()

    np.savetxt(f"{args.outdir}/pdb_entry_nodes.csv.gz", pdb_nodes, delimiter='\t',fmt='%s', header='pdbEntry:ID(pdb-id)',comments='')

    cath_protein_entities = cath_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]
    scop_protein_entities = scop_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]
    interpro_protein_entities = interpro_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]

    protein_entities = pd.concat([cath_protein_entities, scop_protein_entities, interpro_protein_entities]).drop_duplicates(subset = ["protein_entity_id", "protein_entity_ec"])
    protein_entities["ec_list"] = protein_entities["ec_list"].str.replace(",", "|")
    protein_entities.loc[(protein_entities.protein_entity_ec.str.contains("-") == False) & (protein_entities.protein_entity_ec != protein_entities.ec_list), "updatedEC"] = "True"
    protein_entities.loc[(protein_entities.protein_entity_ec.str.contains("-")), "partialEC"] = "True"
    protein_entities["partialEC"] = protein_entities["partialEC"].fillna("False")
    protein_entities["updatedEC"] = protein_entities["updatedEC"].fillna("False")
    protein_entities.rename(columns = {"protein_entity_id": "pdbProteinChain:ID(pdbp-id)", "ec_list" : "ecList:string[]", "protein_entity_ec": "originalEC"}, inplace = True)
    protein_entities.to_csv(f"{args.outdir}/pdb_protein_chain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    protein_ec_rels = protein_entities[["pdbProteinChain:ID(pdbp-id)","ecList:string[]"]].copy()
    protein_ec_rels["ecList:string[]"] = protein_ec_rels["ecList:string[]"].str.split("|")
    protein_ec_rels = protein_ec_rels.explode("ecList:string[]")
    protein_ec_rels = protein_ec_rels.loc[(protein_ec_rels["ecList:string[]"] != "") & (protein_ec_rels["ecList:string[]"].isna() == False)].reset_index()
    protein_ec_rels = protein_ec_rels.drop_duplicates()
    #need to work out why eclist process is returning duplicates
    protein_ec_rels.rename(columns = {"pdbProteinChain:ID(pdbp-id)": ":START_ID(pdbp-id)", "ecList:string[]": ":END_ID(ec-id)"}, inplace = True)
    protein_ec_rels.to_csv(f"{args.outdir}/protein_ec_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domains_nodes = scop_domains[["scop_id", "dm_description"]].drop_duplicates()
    scop_domains_nodes["type"] = "SCOP"
    scop_domains_nodes.rename(columns = {"scop_id": "domain:ID(domain-id)", "dm_description": "name"}, inplace = True)
    #scop_domains_nodes.to_csv(f"{args.outdir}/scop_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_domains_nodes = cath_domains[["cath_domain", "cath_name"]].drop_duplicates()
    cath_domains_nodes["type"] = "CATH"
    cath_domains_nodes.rename(columns = {"cath_domain": "domain:ID(domain-ID)", "cath_name": "name"}, inplace = True)
    #cath_domains_nodes.to_csv(f"{args.outdir}/cath_domains_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    interpro_domain_nodes = interpro_domains[["interpro_accession", "interpro_name", "interpro_type"]].drop_duplicates()
    interpro_domain_nodes["type"] = "InterPro"
    interpro_domain_nodes.rename(columns = {"interpro_accession": "domain:ID(domain-ID)", "interpro_name": "name", "interpro_type": "interproType"}, inplace = True)
    #interpro_domain_nodes.to_csv(f"{args.outdir}/interpro_domain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    domain_nodes = pd.concat([scop_domains_nodes, cath_domains_nodes, interpro_domain_nodes])
    domain_nodes.to_csv(f"{args.outdir}/domain_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_family_nodes = scop_domains[["scop_sunid", "sf_description"]].drop_duplicates()
    scop_family_nodes.rename(columns = {"scop_sunid": "scopFamily:ID(scop-family-id)", "sf_description": "familyDescription"}, inplace = True)
    scop_family_nodes.to_csv(f"{args.outdir}/scop_family_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_superfamily_nodes = scop_domains[["sf_id", "sf_description"]].drop_duplicates()
    scop_superfamily_nodes.rename(columns = {"sf_id": "scopSuperfamily:ID(scop-superfam-id)", "sf_description": "superfamilyDescription"}, inplace = True)
    scop_superfamily_nodes.to_csv(f"{args.outdir}/scop_superfamily_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_class_nodes = scop_domains[["cl_id", "cl_description"]].drop_duplicates()
    scop_class_nodes.rename(columns = {"cl_id": "scopClass:ID(scop-class-id)", "cl_description": "classDescription"}, inplace = True)
    scop_class_nodes.to_csv(f"{args.outdir}/scop_class_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_fold_nodes = scop_domains[["cf_id", "cf_description"]].drop_duplicates()
    scop_fold_nodes.rename(columns = {"cf_id": "scopFold:ID(scop-fold-id)", "cf_description": "foldDescription"}, inplace = True)
    scop_fold_nodes.to_csv(f"{args.outdir}/scop_fold_nodes.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domain_family_rels = scop_domains[["scop_id", "scop_sunid"]].drop_duplicates()
    scop_domain_family_rels.rename(columns = {"scop_id": ":START_ID(domain-id)", "scop_sunid": ":END_ID(scop-family-id)"}, inplace = True)
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
    cath_homology_nodes = cath_domains.cath_homology.unique()

    np.savetxt(f"{args.outdir}/cath_class_nodes.csv.gz", cath_class_nodes, delimiter='\t',fmt='%s', header='cathClass:ID(cath-class-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_architecture_nodes.csv.gz", cath_architecture_nodes, delimiter='\t',fmt='%s', header='cathArchitecture:ID(cath-architecture-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_topology_nodes.csv.gz", cath_topology_nodes, delimiter='\t',fmt='%s', header='cathTopology:ID(cath-topology-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_homology_nodes.csv.gz", cath_homology_nodes, delimiter='\t',fmt='%s', header='cathHomology:ID(cath-homology-ID)', comments='')

    cath_class_architecture_rels = cath_domains[["cath_class", "cath_architecture"]].rename(columns = {"cath_class": ":END_ID(cath-class-ID)", "cath_architecture" : ":START_ID(cath-architecture-ID)"}).drop_duplicates()
    cath_architecture_topology_rels = cath_domains[["cath_architecture", "cath_topology"]].rename(columns = {"cath_architecture": ":END_ID(cath-architecture-ID)", "cath_topology" : ":START_ID(cath-topology-ID)"}).drop_duplicates()
    cath_topology_homology_rels = cath_domains[["cath_topology", "cath_homology"]].rename(columns = {"cath_topology": ":END_ID(cath-topology-ID)", "cath_homology" : ":START_ID(cath-homology-ID)"}).drop_duplicates()
    cath_homology_domain_rels = cath_domains[["cath_homology", "cath_domain"]].rename(columns = {"cath_domain": ":START_ID(domain-ID)", "cath_homology" : ":END_ID(cath-homology-ID)"}).drop_duplicates()

    cath_class_architecture_rels.to_csv(f"{args.outdir}/cath_class_architecture_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_architecture_topology_rels.to_csv(f"{args.outdir}/cath_architecture_topology_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_topology_homology_rels.to_csv(f"{args.outdir}/cath_topology_homology_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    cath_homology_domain_rels.to_csv(f"{args.outdir}/cath_homology_domain_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    bound_molecules = pd.concat([cath_domains.bound_molecule_id, scop_domains.bound_molecule_id, interpro_domains.bound_molecule_id]).drop_duplicates()
    bound_molecules.rename("boundMolecule:ID(bm-id)", inplace = True)
    bound_molecules.to_csv(f"{args.outdir}/bound_molecules.csv.gz", compression = "gzip", sep = "\t", index = False)

    bm_pdb_rels = pd.concat([cath_domains[["bound_molecule_id", "pdb_id"]], scop_domains[["bound_molecule_id", "pdb_id"]], interpro_domains[["bound_molecule_id", "pdb_id"]]]).drop_duplicates()
    bm_pdb_rels.rename(columns = {"bound_molecule_id": ":START_ID(bm-id)", "pdb_id": ":END_ID(pdb-id)"}, inplace = True)
    bm_pdb_rels.to_csv(f"{args.outdir}/bound_molecules_pdb_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    bound_ligand_descriptors = pd.read_pickle(f"{args.bound_ligand_descriptors}")
    bound_sugar_descriptors = pd.read_pickle(f"{args.bound_molecules_sugars_smiles}")
    bound_sugar_descriptors = bound_sugar_descriptors.loc[bound_sugar_descriptors.descriptor.isna() == False]
    bound_sugar_descriptors["ligand_index"] = bound_sugar_descriptors.ligand_index.astype("int") #now that na are removed can coerce to int
    
    cath_bound_ligands = cath_domains.loc[cath_domains.type == "ligand", ["uniqueID", "bound_ligand_id", "name", "type"]].drop_duplicates()
    scop_bound_ligands = scop_domains.loc[scop_domains.type == "ligand",["uniqueID", "bound_ligand_id", "name", "type"]].drop_duplicates()
    interpro_bound_ligands = interpro_domains.loc[interpro_domains.type == "ligand", ["uniqueID","bound_ligand_id", "name", "type"]].drop_duplicates()

    bound_ligands = pd.concat([cath_bound_ligands, scop_bound_ligands, interpro_bound_ligands]).drop_duplicates()
    bound_ligands = bound_ligands.merge(bound_ligand_descriptors, left_on = "bound_ligand_id", right_on = "bl_id", how = "left", indicator = True)
    assert(len(bound_ligands.loc[bound_ligands._merge != "both"]) == 0)
    bound_ligands = bound_ligands.groupby(["uniqueID", "bound_ligand_id", "name", "ligand_entity_description", "descriptor", "ligand_index", "type"]).agg({"ec_list": set}).reset_index()
    bound_ligands["bound_ligand_id"] = bound_ligands["bound_ligand_id"].apply(lambda x: [x]).str.join("|") #convert the string to a list to reflect sugar bound_ligand_id representation (for ligands it is always a single id)
    bound_ligands.rename(columns = {"ligand_entity_description": "entityName", "name" : "entityID", "bound_ligand_id": "componentLigands:string[]", "uniqueID" : "uniqueID:ID(be-id)"}, inplace = True)
    
    cath_bound_sugars = cath_domains.loc[cath_domains.type == "sugar", ["uniqueID", "name", "bound_ligand_id", "type"]].drop_duplicates()
    scop_bound_sugars = scop_domains.loc[scop_domains.type == "sugar", ["uniqueID", "name", "bound_ligand_id", "type"]].drop_duplicates()
    interpro_bound_sugars = interpro_domains.loc[interpro_domains.type == "sugar", ["uniqueID", "name", "bound_ligand_id", "type"]].drop_duplicates()

    bound_sugars = pd.concat([cath_bound_sugars, scop_bound_sugars, interpro_bound_sugars]).drop_duplicates()
    bound_sugars = bound_sugars.groupby(["uniqueID", "name", "type"]).agg({"bound_ligand_id": list}).reset_index()
    bound_sugars["bound_ligand_id"] = bound_sugars["bound_ligand_id"].str.join("|")
    bound_sugars = bound_sugars.merge(bound_sugar_descriptors[["uniqueID", "descriptor", "ec_list", "ligand_index"]], on = "uniqueID", how = "left", indicator = True)
    assert(len(bound_sugars.loc[bound_sugars._merge != "both"]) == 0)
    bound_sugars = bound_sugars.groupby(["uniqueID", "name", "bound_ligand_id", "descriptor", "ligand_index", "type"]).agg({"ec_list": set}).reset_index()
    bound_sugars["entityID"] = bound_sugars["name"]

    bound_sugars.rename(columns = {"uniqueID": "uniqueID:ID(be-id)", "name": "entityName", "bound_ligand_id": "componentLigands:string[]"}, inplace = True)
    
    bound_entities = pd.concat([bound_sugars, bound_ligands])
    bound_entities["ec_list"] = bound_entities.ec_list.str.join("|")
    bound_entities.rename(columns = {"ec_list": "ecList:string[]"}, inplace = True)
    bound_entities.to_csv(f"{args.outdir}/bound_entities.csv.gz", compression = "gzip", sep = "\t", index = False)

    bm_be_rels = pd.concat(
        [cath_domains[["bound_molecule_id", "uniqueID"]],
         scop_domains[["bound_molecule_id", "uniqueID"]], 
         interpro_domains[["bound_molecule_id", "uniqueID"]]
        ]).drop_duplicates().rename(columns = {"uniqueID": ":END_ID(be-id)", "bound_molecule_id": ":START_ID(bm-id)"})

    bm_be_rels.to_csv(f"{args.outdir}/bm_be_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    parity_calcs = pd.read_pickle(f"{args.parity_calcs}")
    parity_calcs = parity_calcs.explode("ec")
    parity_calcs_filtered = parity_calcs.loc[parity_calcs.score >= 0.25]

    parity_rels = bound_entities[["componentLigands:string[]", "uniqueID:ID(be-id)", "ligand_index", "ecList:string[]"]].copy()
    parity_rels["ecList:string[]"] = parity_rels["ecList:string[]"].str.split("|")
    parity_rels = parity_rels.explode("ecList:string[]")
    parity_rels = parity_rels.merge(parity_calcs_filtered, left_on = ["ligand_index", "ecList:string[]"], right_on = ["pdb_ligand", "ec"], how = "inner")
    parity_rels = parity_rels.merge(biological_ligands[["canonical_smiles", "uniqueID"]].drop_duplicates(), left_on = "compound", right_on = "canonical_smiles", how = "left", indicator = True)
    assert(len(parity_rels.loc[parity_rels._merge != "both"]) == 0)
    parity_rels.drop(columns = ["_merge"], inplace = True)
    parity_rels = parity_rels[["uniqueID:ID(be-id)", "uniqueID", "pdbl_subparity", "score", "ec"]].rename(columns = {"uniqueID:ID(be-id)": ":START_ID(be-id)", "uniqueID": ":END_ID(bio-id)"})
    parity_rels = parity_rels.groupby([col for col in parity_rels.columns if col != "ec"]).agg({"ec": list}).reset_index()
    parity_rels["ec"] = parity_rels["ec"].str.join("|")
    parity_rels.rename(columns = {"score": "parityScore:float", "pdbl_subparity": "subParityScore:float", "ec": "ecList:string[]"}, inplace = True)
    parity_rels.to_csv(f"{args.outdir}/bound_entity_parity_score_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_domain_ligand_interactions = cath_domains[["cath_domain", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "uniqueID"]].drop_duplicates()
    cath_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "cath_domain": ":START_ID(domain-ID)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc" , "domain_ownership" : "interactionMode"}, inplace = True)
    cath_domain_ligand_interactions.to_csv(f"{args.outdir}/cath_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_domain_ligand_interactions = scop_domains[["scop_id", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "uniqueID"]].drop_duplicates()
    scop_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "scop_id": ":START_ID(domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_ownership" : "interactionMode"}, inplace = True)
    scop_domain_ligand_interactions.to_csv(f"{args.outdir}/scop_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)

    interpro_domain_ligand_interactions = interpro_domains[["interpro_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "uniqueID"]].drop_duplicates()
    interpro_domain_ligand_interactions.rename(columns = {"uniqueID": ":END_ID(be-id)", "interpro_accession": ":START_ID(domain-ID)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_ownership" : "interactionMode"}, inplace = True)
    interpro_domain_ligand_interactions.to_csv(f"{args.outdir}/interpro_domain_ligand_interactions.csv.gz", compression = "gzip", sep = "\t", index = False)

    scop_pdb_protein_rels = scop_domains[["pdb_id", "protein_entity_id"]].drop_duplicates()
    cath_pdb_protein_rels = cath_domains[["pdb_id", "protein_entity_id"]].drop_duplicates()
    interpro_pdb_protein_rels = interpro_domains[["pdb_id", "protein_entity_id"]].drop_duplicates()

    pdb_protein_rels = pd.concat([scop_pdb_protein_rels, cath_pdb_protein_rels, interpro_pdb_protein_rels]).drop_duplicates()
    pdb_protein_rels.rename(columns = {"protein_entity_id": ":START_ID(pdbp-id)", "pdb_id": ":END_ID(pdb-id)"}, inplace = True)
    pdb_protein_rels.to_csv(f"{args.outdir}/pdb_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_protein_rels = cath_domains[["cath_domain", "protein_entity_id"]].drop_duplicates().rename(columns = {"cath_domain": ":START_ID(domain-ID)", "protein_entity_id":":END_ID(pdbp-id)"})
    scop_protein_rels = scop_domains[["scop_id", "protein_entity_id"]].drop_duplicates().rename(columns = {"scop_id": ":START_ID(domain-id)", "protein_entity_id":":END_ID(pdbp-id)"})
    interpro_protein_rels = interpro_domains[["interpro_accession", "protein_entity_id"]].drop_duplicates().rename(columns = {"interpro_accession": ":START_ID(domain-ID)", "protein_entity_id":":END_ID(pdbp-id)"})
    
    cath_protein_rels.to_csv(f"{args.outdir}/cath_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    scop_protein_rels.to_csv(f"{args.outdir}/scop_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)
    interpro_protein_rels.to_csv(f"{args.outdir}/interpro_protein_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    procoggraph_node = pd.DataFrame({"procoggraph:ID(procoggraph-id)": ["procoggraph"],
                                     "name": ["ProCogGraph"],
                                     "description": ["procoggraph"],
                                     "date_created": ["2023"],
                                     "date_updated": ["2023"],
                                     "database_version": ["0.1"],
                                     "biological_ligands_version": ["0.1"],
                                     "pdbe_graph_version": ["0.1"],
                                     "pdbe_graph_scripts_version": ["0.1"],
                                     "pdbe_graph_data_version": ["0.1"],
                                     "input_params": ["-"],})
    
    procoggraph_node.to_csv(f"{args.outdir}/procoggraph_node.csv.gz", compression = "gzip", sep = "\t", index = False)

if __name__ == "__main__":
    main()