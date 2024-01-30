#!/usr/bin/env python

import pandas as pd
from rdkit import Chem
import numpy as np
import re
from Bio.ExPASy import Enzyme as EEnzyme
from pathlib import Path
import argparse
from utils import get_terminal_record, process_ec_records

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

    return ec_records_df_grouped, ec_class_descriptions, ec_subclass_descriptions, ec_subsubclass_descriptions

#make an argument parser with the following arguments: enzdat, enzclass, biological_ligands, domain_ownership, bound_ligand_descriptors, bound_molecules_sugars_ec, parity_calcs

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
    parser.add_argument('--domain_ownership', metavar='domain_ownership', type=str,
                        help = "path to domain ownership dataframe")
    parser.add_argument('--bound_ligand_descriptors', metavar='bound_ligand_descriptors', type=str,
                        help = "path to bound ligand descriptors dataframe")
    parser.add_argument('--bound_molecules_sugars_ec', metavar='bound_molecules_sugars_ec', type=str,
                        help = "path to bound molecules sugars ec dataframe")

    args = parser.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    #check if all files to be produced already exist in dir, exit if they do
    #if 

    ec_records_df_grouped = process_ec_records(args.enzyme_dat_file , args.enzyme_class_file)


    ec_id_nodes = ec_records_df_grouped[["TRANSFER", "DE"]].rename(columns = {"TRANSFER" : "ecID:ID(ec-id)", "DE" : "description"}).drop_duplicates()
    ec_nodes_class = ec_records_df_grouped[["class", "class_description"]].rename(columns = {"class": "ecClass:ID(class-id)", "class_description": "description"}).drop_duplicates()
    ec_nodes_subclass = ec_records_df_grouped[["subclass", "subclass_description"]].rename(columns = {"subclass": "ecSubClass:ID(subclass-id)", "subclass_description": "description"}).drop_duplicates()
    ec_nodes_subsubclass = ec_records_df_grouped[["subsubclass", "subsubclass_description"]].rename(columns = {"subsubclass": "ecSubSubClass:ID(subsubclass-id)", "subsubclass_description": "description"}).drop_duplicates()

    ec_id_nodes.to_csv(f"{args.outdir}/ec_id_nodes.csv", sep = "\t", index = False)
    ec_nodes_class.to_csv(f"{args.outdir}/ec_nodes_class.csv", sep = "\t", index = False)
    ec_nodes_subclass.to_csv(f"{args.outdir}/ec_nodes_subclass.csv", sep = "\t", index = False)
    ec_nodes_subsubclass.to_csv(f"{args.outdir}/ec_nodes_subsubclass.csv", sep = "\t", index = False)

    ec_subsubclass_id_rel = ec_records_df_grouped[["TRANSFER", "subsubclass"]].drop_duplicates()
    ec_subsubclass_id_rel.rename(columns = {"TRANSFER" : ":START_ID(ec-id)", "subsubclass": ":END_ID(subsubclass-id)"}, inplace = True)
    ec_subclass_subsubclass_rel = ec_records_df_grouped[["subclass", "subsubclass"]].drop_duplicates()
    ec_subclass_subsubclass_rel.rename(columns = {"subsubclass" : ":START_ID(subsubclass-id)", "subclass": ":END_ID(subclass-id)"}, inplace = True)
    ec_class_subclass_rel = ec_records_df_grouped[["class", "subclass"]].drop_duplicates()
    ec_class_subclass_rel.rename(columns = {"subclass" : ":START_ID(subclass-id)", "class": ":END_ID(class-id)"}, inplace = True)

    ec_class_subclass_rel.to_csv(f"{args.outdir}/ec_class_subclass_rel.csv", sep = "\t", index = False)
    ec_subclass_subsubclass_rel.to_csv(f"{args.outdir}/ec_subclass_subsubclass_rel.csv", sep = "\t", index = False)
    ec_subsubclass_id_rel.to_csv(f"{args.outdir}/ec_subsubclass_id_rel.csv", sep = "\t", index = False)

    #It is questionable whether 204 properties is a good idea for a node. how can this be refactored out into nodes in the graph?
    biological_ligands = pd.read_pickle(args.biological_ligands)

    biological_ligands_nodes = biological_ligands[[col for col in biological_ligands.columns if col not in ["entry", "_merge"]]].copy()
    biological_ligands_nodes.drop_duplicates(subset = ["canonical_smiles", "uniqueID"], inplace = True)
    biological_ligands_nodes = biological_ligands_nodes[["canonical_smiles", "uniqueID", "ligand_db"]]
    biological_ligands_nodes.rename(columns = {"uniqueID": "uniqueID:ID(bio-id)", "canonical_smiles": "canonicalSMILES", "ligand_db" : "ligandDB:string[]"}, inplace = True)

    biological_ligands_nodes.to_csv(f"{args.outdir}/biological_ligand_nodes.csv",sep = "\t", index = False)

    biological_ligands_ec = biological_ligands[["uniqueID", "entry"]].drop_duplicates()
    biological_ligands_ec.rename(columns = {"uniqueID": ":START_ID(bio-id)", "entry": ":END_ID(ec-id)"}, inplace = True)
    biological_ligands_ec.to_csv(f"{args.outdir}/biological_ligands_ec.csv", sep = "\t", index = False)

    cath_bl_domains = pd.read_csv("../domain_ownership/cath_bl_domain_ownership.csv", na_values = ["NaN", "None"], keep_default_na = False)
    scop_bl_domains = pd.read_csv("../domain_ownership/scop_bl_domain_ownership.csv", na_values = ["NaN", "None"], keep_default_na = False)
    interpro_bl_domains = pd.read_csv("../domain_ownership/interpro_bl_domain_ownership.csv", na_values = ["NaN", "None"], keep_default_na = False)

    cath_bl_domains["bound_ligand_unique_id"] = cath_bl_domains["bound_molecule_id"] + "_" + cath_bl_domains["bound_ligand_id"]
    scop_bl_domains["bound_ligand_unique_id"] = scop_bl_domains["bound_molecule_id"] + "_" + scop_bl_domains["bound_ligand_id"]
    interpro_bl_domains["bound_ligand_unique_id"] = interpro_bl_domains["bound_molecule_id"] + "_" + interpro_bl_domains["bound_ligand_id"]

    cath_sugar_domains = pd.read_csv("../domain_ownership/cath_sugar_domain_ownership.csv", na_values = ["NaN", "None"], keep_default_na = False)
    scop_sugar_domains = pd.read_csv("../domain_ownership/scop_sugar_domain_ownership.csv", na_values = ["NaN", "None"], keep_default_na = False)
    interpro_sugar_domains = pd.read_csv("../domain_ownership/interpro_sugar_domain_ownership.csv", na_values = ["NaN", "None"], keep_default_na = False)

    pdb_nodes = pd.concat([cath_bl_domains.pdb_id, scop_bl_domains.pdb_id, interpro_bl_domains.pdb_id,
                        cath_sugar_domains.pdb_id, scop_sugar_domains.pdb_id, interpro_sugar_domains.pdb_id]).unique()

    np.savetxt(f"{args.outdir}/pdb_entry_nodes.csv", pdb_nodes, delimiter='\t',fmt='%s', header='pdbEntry:ID(pdb-id)',comments='')

    cath_bl_protein_entities = cath_bl_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]
    scop_bl_protein_entities = scop_bl_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]
    interpro_bl_protein_entities = interpro_bl_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]

    cath_sugar_protein_entities = cath_sugar_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]
    scop_sugar_protein_entities = scop_sugar_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]
    interpro_sugar_protein_entities = interpro_sugar_domains[["protein_entity_id", "protein_entity_ec", "ec_list"]]

    protein_entities = pd.concat([cath_bl_protein_entities, scop_bl_protein_entities, interpro_bl_protein_entities,
                                cath_sugar_protein_entities, scop_sugar_protein_entities, interpro_sugar_protein_entities]).drop_duplicates()
    protein_entities["ec_list"] = protein_entities["ec_list"].str.replace(",", "|")
    protein_entities.loc[(protein_entities.protein_entity_ec.str.contains("-") == False) & (protein_entities.protein_entity_ec != protein_entities.ec_list), "updatedEC"] = "True"
    protein_entities.loc[(protein_entities.protein_entity_ec.str.contains("-")), "partialEC"] = "True"
    protein_entities["partialEC"] = protein_entities["partialEC"].fillna("False")
    protein_entities["updatedEC"] = protein_entities["updatedEC"].fillna("False")
    protein_entities.rename(columns = {"protein_entity_id": "pdbProteinChain:ID(pdbp-id)", "ec_list" : "ecList:string[]", "protein_entity_ec": "originalEC"}, inplace = True)
    protein_entities.to_csv(f"{args.outdir}/pdb_protein_chain_nodes.csv", sep = "\t", index = False)

    protein_ec_rels = protein_entities[["pdbProteinChain:ID(pdbp-id)","ecList:string[]"]].copy()
    protein_ec_rels["ecList:string[]"] = protein_ec_rels["ecList:string[]"].str.split("|")
    protein_ec_rels = protein_ec_rels.explode("ecList:string[]")
    protein_ec_rels = protein_ec_rels.loc[(protein_ec_rels["ecList:string[]"] != "") & (protein_ec_rels["ecList:string[]"].isna() == False)].reset_index()
    protein_ec_rels = protein_ec_rels.drop_duplicates()
    #need to work out why eclist process is returning duplicates
    protein_ec_rels.rename(columns = {"pdbProteinChain:ID(pdbp-id)": ":START_ID(pdbp-id)", "ecList:string[]": ":END_ID(ec-id)"}, inplace = True)
    protein_ec_rels.to_csv(f"{args.outdir}/protein_ec_rels.csv", sep = "\t", index = False)

    scop_bl_domains_nodes = scop_bl_domains[["scop_id", "dm_description"]].drop_duplicates()
    scop_sugar_domains_nodes = scop_sugar_domains[["scop_id", "dm_description"]].drop_duplicates()
    scop_domains_nodes = pd.concat([scop_bl_domains_nodes, scop_sugar_domains_nodes]).drop_duplicates()
    scop_domains_nodes.rename(columns = {"scop_id": "scopDomain:ID(scop-domain-id)", "dm_description": "domainDescription"}, inplace = True)
    scop_domains_nodes.to_csv(f"{args.outdir}/scop_domains_nodes.csv", sep = "\t", index = False)

    scop_bl_family_nodes = scop_bl_domains[["scop_sunid", "sf_description"]].drop_duplicates()
    scop_sugar_family_nodes = scop_sugar_domains[["scop_sunid", "sf_description"]].drop_duplicates()
    scop_family_nodes = pd.concat([scop_bl_family_nodes,scop_sugar_family_nodes]).drop_duplicates()
    scop_family_nodes.rename(columns = {"scop_sunid": "scopFamily:ID(scop-family-id)", "sf_description": "familyDescription"}, inplace = True)
    scop_family_nodes.to_csv(f"{args.outdir}/scop_family_nodes.csv", sep = "\t", index = False)

    scop_bl_superfamily_nodes = scop_bl_domains[["sf_id", "sf_description"]].drop_duplicates()
    scop_sugar_superfamily_nodes = scop_sugar_domains[["sf_id", "sf_description"]].drop_duplicates()
    scop_superfamily_nodes = pd.concat([scop_bl_superfamily_nodes, scop_sugar_superfamily_nodes]).drop_duplicates()
    scop_superfamily_nodes.rename(columns = {"sf_id": "scopSuperfamily:ID(scop-superfam-id)", "sf_description": "superfamilyDescription"}, inplace = True)
    scop_superfamily_nodes.to_csv(f"{args.outdir}/scop_superfamily_nodes.csv", sep = "\t", index = False)

    scop_bl_class_nodes = scop_bl_domains[["cl_id", "cl_description"]].drop_duplicates()
    scop_sugar_class_nodes = scop_sugar_domains[["cl_id", "cl_description"]].drop_duplicates()
    scop_class_nodes = pd.concat([scop_bl_class_nodes, scop_sugar_class_nodes]).drop_duplicates()
    scop_class_nodes.rename(columns = {"cl_id": "scopClass:ID(scop-class-id)", "cl_description": "classDescription"}, inplace = True)
    scop_class_nodes.to_csv(f"{args.outdir}/scop_class_nodes.csv", sep = "\t", index = False)

    scop_bl_fold_nodes = scop_bl_domains[["cf_id", "cf_description"]].drop_duplicates()
    scop_sugar_fold_nodes = scop_sugar_domains[["cf_id", "cf_description"]].drop_duplicates()
    scop_fold_nodes = pd.concat([scop_bl_fold_nodes,scop_sugar_fold_nodes]).drop_duplicates()
    scop_fold_nodes.rename(columns = {"cf_id": "scopFold:ID(scop-fold-id)", "cf_description": "foldDescription"}, inplace = True)
    scop_fold_nodes.to_csv(f"{args.outdir}/scop_fold_nodes.csv", sep = "\t", index = False)

    scop_bl_domain_family_rels = scop_bl_domains[["scop_id", "scop_sunid"]].drop_duplicates()
    scop_sugar_domain_family_rels = scop_sugar_domains[["scop_id", "scop_sunid"]].drop_duplicates()
    scop_domain_family_rels = pd.concat([scop_bl_domain_family_rels,scop_sugar_domain_family_rels]).drop_duplicates()
    scop_domain_family_rels.rename(columns = {"scop_id": ":START_ID(scop-domain-id)", "scop_sunid": ":END_ID(scop-family-id)"}, inplace = True)
    scop_domain_family_rels.to_csv(f"{args.outdir}/scop_domain_family_rels.csv", sep = "\t", index = False)

    scop_bl_family_superfamily_rels = scop_bl_domains[["scop_sunid", "sf_id"]].drop_duplicates()
    scop_sugar_family_superfamily_rels = scop_sugar_domains[["scop_sunid", "sf_id"]].drop_duplicates()
    scop_family_superfamily_rels = pd.concat([scop_bl_family_superfamily_rels, scop_sugar_family_superfamily_rels]).drop_duplicates()
    scop_family_superfamily_rels.rename(columns = {"scop_sunid": ":START_ID(scop-family-id)", "sf_id": ":END_ID(scop-superfam-id)"}, inplace = True)
    scop_family_superfamily_rels.to_csv(f"{args.outdir}/scop_family_superfam_rels.csv", sep = "\t", index = False)

    scop_bl_superfamily_fold_rels = scop_bl_domains[["sf_id", "cf_id"]].drop_duplicates()
    scop_sugar_superfamily_fold_rels = scop_sugar_domains[["sf_id", "cf_id"]].drop_duplicates()
    scop_superfamily_fold_rels = pd.concat([scop_bl_superfamily_fold_rels,scop_sugar_superfamily_fold_rels]).drop_duplicates()
    scop_superfamily_fold_rels.rename(columns = {"sf_id": ":START_ID(scop-superfam-id)", "cf_id": ":END_ID(scop-fold-id)"}, inplace = True)
    scop_superfamily_fold_rels.to_csv(f"{args.outdir}/scop_superfam_fold_rels.csv", sep = "\t", index = False)

    scop_bl_fold_class_rels = scop_bl_domains[["cf_id", "cl_id"]].drop_duplicates()
    scop_sugar_fold_class_rels = scop_sugar_domains[["cf_id", "cl_id"]].drop_duplicates()
    scop_fold_class_rels = pd.concat([scop_bl_fold_class_rels, scop_sugar_fold_class_rels]).drop_duplicates()
    scop_fold_class_rels.rename(columns = {"cf_id": ":START_ID(scop-fold-id)", "cl_id": ":END_ID(scop-class-id)"}, inplace = True)
    scop_fold_class_rels.to_csv(f"{args.outdir}/scop_fold_class_rels.csv", sep = "\t", index = False)

    cath_bl_domains_nodes = cath_bl_domains[["cath_domain", "cath_name"]].drop_duplicates()
    cath_sugar_domains_nodes = cath_sugar_domains[["cath_domain", "cath_name"]].drop_duplicates()
    cath_domains_nodes = pd.concat([cath_bl_domains_nodes,cath_sugar_domains_nodes]).drop_duplicates()
    cath_domains_nodes.rename(columns = {"cath_domain": "cathDomain:ID(cath-domain-ID)", "cath_name": "cathName"}, inplace = True)
    cath_domains_nodes.to_csv(f"{args.outdir}/cath_domains_nodes.csv", sep = "\t", index = False)

    cath_class_nodes = pd.concat([cath_bl_domains.cath_class, cath_sugar_domains.cath_class]).unique()
    cath_architecture_nodes = pd.concat([cath_bl_domains.cath_architecture, cath_sugar_domains.cath_architecture]).unique()
    cath_topology_nodes = pd.concat([cath_bl_domains.cath_topology,cath_sugar_domains.cath_topology]).unique()
    cath_homology_nodes = pd.concat([cath_bl_domains.cath_homology,cath_sugar_domains.cath_homology]).unique()

    np.savetxt(f"{args.outdir}/cath_class_nodes.csv", cath_class_nodes, delimiter='\t',fmt='%s', header='cathClass:ID(cath-class-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_architecture_nodes.csv", cath_architecture_nodes, delimiter='\t',fmt='%s', header='cathArchitecture:ID(cath-architecture-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_topology_nodes.csv", cath_topology_nodes, delimiter='\t',fmt='%s', header='cathTopology:ID(cath-topology-ID)', comments='')
    np.savetxt(f"{args.outdir}/cath_homology_nodes.csv", cath_homology_nodes, delimiter='\t',fmt='%s', header='cathHomology:ID(cath-homology-ID)', comments='')

    cath_bl_class_architecture_rels = cath_bl_domains[["cath_class", "cath_architecture"]].rename(columns = {"cath_class": ":END_ID(cath-class-ID)", "cath_architecture" : ":START_ID(cath-architecture-ID)"}).drop_duplicates()
    cath_sugar_class_architecture_rels = cath_sugar_domains[["cath_class", "cath_architecture"]].rename(columns = {"cath_class": ":END_ID(cath-class-ID)", "cath_architecture" : ":START_ID(cath-architecture-ID)"}).drop_duplicates()
    cath_class_architecture_rels = pd.concat([cath_bl_class_architecture_rels, cath_sugar_class_architecture_rels]).drop_duplicates()

    cath_bl_architecture_topology_rels = cath_bl_domains[["cath_architecture", "cath_topology"]].rename(columns = {"cath_architecture": ":END_ID(cath-architecture-ID)", "cath_topology" : ":START_ID(cath-topology-ID)"}).drop_duplicates()
    cath_sugar_architecture_topology_rels = cath_sugar_domains[["cath_architecture", "cath_topology"]].rename(columns = {"cath_architecture": ":END_ID(cath-architecture-ID)", "cath_topology" : ":START_ID(cath-topology-ID)"}).drop_duplicates()
    cath_architecture_topology_rels = pd.concat([cath_bl_architecture_topology_rels,cath_sugar_architecture_topology_rels]).drop_duplicates()

    cath_bl_topology_homology_rels = cath_bl_domains[["cath_topology", "cath_homology"]].rename(columns = {"cath_topology": ":END_ID(cath-topology-ID)", "cath_homology" : ":START_ID(cath-homology-ID)"}).drop_duplicates()
    cath_sugar_topology_homology_rels = cath_sugar_domains[["cath_topology", "cath_homology"]].rename(columns = {"cath_topology": ":END_ID(cath-topology-ID)", "cath_homology" : ":START_ID(cath-homology-ID)"}).drop_duplicates()
    cath_topology_homology_rels = pd.concat([cath_bl_topology_homology_rels, cath_sugar_topology_homology_rels]).drop_duplicates()

    cath_bl_homology_domain_rels = cath_bl_domains[["cath_homology", "cath_domain"]].rename(columns = {"cath_domain": ":START_ID(cath-domain-ID)", "cath_homology" : ":END_ID(cath-homology-ID)"}).drop_duplicates()
    cath_sugar_homology_domain_rels = cath_sugar_domains[["cath_homology", "cath_domain"]].rename(columns = {"cath_domain": ":START_ID(cath-domain-ID)", "cath_homology" : ":END_ID(cath-homology-ID)"}).drop_duplicates()
    cath_homology_domain_rels = pd.concat([cath_bl_homology_domain_rels,cath_sugar_homology_domain_rels]).drop_duplicates()

    cath_class_architecture_rels.to_csv(f"{args.outdir}/cath_class_architecture_rels.csv", sep = "\t", index = False)
    cath_architecture_topology_rels.to_csv(f"{args.outdir}/cath_architecture_topology_rels.csv", sep = "\t", index = False)
    cath_topology_homology_rels.to_csv(f"{args.outdir}/cath_topology_homology_rels.csv", sep = "\t", index = False)
    cath_homology_domain_rels.to_csv(f"{args.outdir}/cath_homology_domain_rels.csv", sep = "\t", index = False)

    interpro_bl_domain_nodes = interpro_bl_domains[["interpro_accession", "interpro_name", "interpro_type"]].drop_duplicates()
    interpro_sugar_domain_nodes = interpro_sugar_domains[["interpro_accession", "interpro_name", "interpro_type"]].drop_duplicates()
    interpro_domain_nodes = pd.concat([interpro_bl_domain_nodes,interpro_sugar_domain_nodes]).drop_duplicates()
    interpro_domain_nodes.rename(columns = {"interpro_accession": "interpro_domain:ID(ipr-domain-ID)", "interpro_name": "interproName", "interpro_type": "interproType"}, inplace = True)
    interpro_domain_nodes.to_csv(f"{args.outdir}/interpro_domain_nodes.csv", sep = "\t", index = False)

    bound_ligand_descriptors = pd.read_pickle("../pdbe_graph_files/bound_ligand_descriptors.pkl")
    bound_sugar_descriptors = pd.read_pickle("../pdbe_graph_files/bound_molecules_sugars_ec.pkl")
    bound_sugar_descriptors = bound_sugar_descriptors.loc[bound_sugar_descriptors.descriptor.isna() == False]
    bound_sugar_descriptors["ligand_index"] = bound_sugar_descriptors.ligand_index.astype("int") #now that na are removed can coerce to int
    #should try and make sure this is created on df creation
    bound_sugar_descriptors["sugar_id"] = bound_sugar_descriptors["bound_molecule_id"] + "_se" + bound_sugar_descriptors["ligand_entity_id_numerical"].astype("str")

    bound_molecules = pd.concat([
        cath_bl_domains.bound_molecule_id, cath_sugar_domains.bound_molecule_id,
        scop_bl_domains.bound_molecule_id, scop_sugar_domains.bound_molecule_id,
        interpro_bl_domains.bound_molecule_id, interpro_sugar_domains.bound_molecule_id]).drop_duplicates()
    bound_molecules.rename("boundMolecule:ID(bm-id)", inplace = True)
    bound_molecules.to_csv(f"{args.outdir}/bound_molecules.csv", sep = "\t", index = False)

    #need to do boundmolecules to pdb entry mapping.

    bm_pdb_rels = pd.concat([cath_bl_domains[["bound_molecule_id", "pdb_id"]], cath_sugar_domains[["bound_molecule_id", "pdb_id"]],
                            scop_bl_domains[["bound_molecule_id", "pdb_id"]], scop_sugar_domains[["bound_molecule_id", "pdb_id"]],
                            interpro_bl_domains[["bound_molecule_id", "pdb_id"]], interpro_sugar_domains[["bound_molecule_id", "pdb_id"]]]).drop_duplicates()
    bm_pdb_rels.rename(columns = {"bound_molecule_id": ":START_ID(bm-id)", "pdb_id": ":END_ID(pdb-id)"}, inplace = True)
    bm_pdb_rels.to_csv(f"{args.outdir}/bound_molecules_pdb_rels.csv", sep = "\t", index = False)

    cath_bound_ligands = cath_bl_domains[["bound_ligand_unique_id", "bound_ligand_id", "bound_ligand_name"]].drop_duplicates()
    scop_bound_ligands = scop_bl_domains[["bound_ligand_unique_id", "bound_ligand_id", "bound_ligand_name"]].drop_duplicates()
    interpro_bound_ligands = interpro_bl_domains[["bound_ligand_unique_id","bound_ligand_id", "bound_ligand_name"]].drop_duplicates()

    bound_ligands = pd.concat([cath_bound_ligands, scop_bound_ligands, interpro_bound_ligands]).drop_duplicates()

    bound_ligands = bound_ligands.merge(bound_ligand_descriptors, left_on = "bound_ligand_id", right_on = "bl_id", how = "left", indicator = True)
    assert(len(bound_ligands.loc[bound_ligands._merge != "both"]) == 0)
    bound_ligands = bound_ligands.groupby(["bound_ligand_unique_id", "bound_ligand_id", "bound_ligand_name", "ligand_entity_description", "descriptor", "ligand_index"]).agg({"ec_list": set}).reset_index()

    cath_bound_sugars = cath_sugar_domains[["sugar_id", "ligand_entity_description", "bound_ligand_id"]].drop_duplicates()
    scop_bound_sugars = scop_sugar_domains[["sugar_id", "ligand_entity_description", "bound_ligand_id"]].drop_duplicates()
    interpro_bound_sugars = interpro_sugar_domains[["sugar_id", "ligand_entity_description", "bound_ligand_id"]].drop_duplicates()

    bound_sugars = pd.concat([cath_bound_sugars, scop_bound_sugars, interpro_bound_sugars]).drop_duplicates()
    bound_sugars = bound_sugars.groupby(["sugar_id", "ligand_entity_description"]).agg({"bound_ligand_id": list}).reset_index()
    bound_sugars["bound_ligand_id"] = bound_sugars["bound_ligand_id"].str.join("|")
    bound_sugars = bound_sugars.merge(bound_sugar_descriptors[["sugar_id", "descriptor", "ec_list", "ligand_index"]], on = "sugar_id", how = "left", indicator = True)
    assert(len(bound_sugars.loc[bound_sugars._merge != "both"]) == 0)
    bound_sugars = bound_sugars.groupby(["sugar_id", "ligand_entity_description", "bound_ligand_id", "descriptor", "ligand_index"]).agg({"ec_list": set}).reset_index()
    bound_sugars["entityID"] = bound_sugars["ligand_entity_description"]

    bound_ligands["type"] = "ligand"
    bound_sugars["type"] = "sugar"

    #bound_ligands["uniqueID:ID(be-id)"] = bound_ligands["bound_ligand_unique_id"]
    bound_ligands["bound_ligand_id"] = bound_ligands["bound_ligand_id"].apply(lambda x: [x]).str.join("|")

    bound_sugars.rename(columns = {"sugar_id": "uniqueID:ID(be-id)", "ligand_entity_description": "entityName", "bound_ligand_id": "componentLigands:string[]"}, inplace = True)
    bound_ligands.rename(columns = {"ligand_entity_description": "entityName", "bound_ligand_name" : "entityID", "bound_ligand_id": "componentLigands:string[]", "bound_ligand_unique_id" : "uniqueID:ID(be-id)"}, inplace = True)

    bound_entities = pd.concat([bound_sugars, bound_ligands])
    bound_entities["ec_list"] = bound_entities.ec_list.str.join("|")
    bound_entities.rename(columns = {"ec_list": "ecList:string[]"}, inplace = True)
    bound_entities.to_csv(f"{args.outdir}/bound_entities.csv", sep = "\t", index = False)

    bm_bl_rels = pd.concat([cath_bl_domains[["bound_molecule_id", "bound_ligand_unique_id"]],
                            scop_bl_domains[["bound_molecule_id", "bound_ligand_unique_id"]], 
                            interpro_bl_domains[["bound_molecule_id", "bound_ligand_unique_id"]]]).drop_duplicates().rename(columns = {"bound_ligand_unique_id": ":END_ID(be-id)", "bound_molecule_id": ":START_ID(bm-id)"})
    bm_sug_rels = pd.concat([cath_sugar_domains[["bound_molecule_id", "sugar_id"]],
                            scop_sugar_domains[["bound_molecule_id", "sugar_id"]],
                            interpro_sugar_domains[["bound_molecule_id", "sugar_id"]]]).drop_duplicates().rename(columns = {"sugar_id": ":END_ID(be-id)", "bound_molecule_id": ":START_ID(bm-id)"})

    bm_be_rels = pd.concat([bm_bl_rels, bm_sug_rels]).drop_duplicates()

    bm_be_rels.to_csv("bm_be_rels.csv", sep = "\t", index = False)

    parity_calcs = pd.read_pickle("../parity_calcs/bound_entities_parity_2/all_parity_calcs.pkl")
    parity_calcs = parity_calcs.explode("ec")
    parity_calcs_filtered = parity_calcs.loc[parity_calcs.score >= 0.25]

    parity_bl = bound_ligands[["componentLigands:string[]", "uniqueID:ID(be-id)", "ligand_index", "ec_list"]].explode("ec_list")
    parity_bl = parity_bl.merge(parity_calcs_filtered, left_on = ["ligand_index", "ec_list"], right_on = ["pdb_ligand", "ec"], how = "inner")
    parity_bl = parity_bl.merge(biological_ligands[["canonical_smiles", "uniqueID"]].drop_duplicates(), left_on = "compound", right_on = "canonical_smiles", how = "left", indicator = True)
    assert(len(parity_bl.loc[parity_bl._merge != "both"]) == 0)
    parity_bl.drop(columns = ["_merge"], inplace = True)
    parity_bl = parity_bl[["uniqueID:ID(be-id)", "uniqueID", "pdbl_subparity", "score", "ec"]].rename(columns = {"uniqueID:ID(be-id)": ":START_ID(be-id)", "uniqueID": ":END_ID(bio-id)"})
    parity_bl = parity_bl.groupby([col for col in parity_bl.columns if col != "ec"]).agg({"ec": list}).reset_index()
    parity_bl["ec"] = parity_bl["ec"].str.join("|")
    parity_bl.rename(columns = {"score": "parityScore:float", "pdbl_subparity": "subParityScore:float", "ec": "ecID:string[]"}, inplace = True)
    #now sugars! §=

    parity_sugars = bound_sugars[["uniqueID:ID(be-id)", "componentLigands:string[]", "ligand_index", "ec_list"]].explode("ec_list")
    parity_sugars = parity_sugars.merge(parity_calcs_filtered, left_on = ["ligand_index", "ec_list"], right_on = ["pdb_ligand", "ec"], how = "inner")
    parity_sugars = parity_sugars.merge(biological_ligands[["canonical_smiles", "uniqueID"]].drop_duplicates(), left_on = "compound", right_on = "canonical_smiles", how = "left", indicator = True)
    assert(len(parity_sugars.loc[parity_sugars._merge != "both"]) == 0)
    parity_sugars.drop(columns = ["_merge"], inplace = True)
    parity_sugars = parity_sugars[["uniqueID:ID(be-id)", "uniqueID", "pdbl_subparity", "score", "ec"]].rename(columns = {"uniqueID:ID(be-id)": ":START_ID(be-id)", "uniqueID": ":END_ID(bio-id)"})
    parity_sugars = parity_sugars.groupby([col for col in parity_sugars.columns if col != "ec"]).agg({"ec": list}).reset_index()
    parity_sugars["ec"] = parity_sugars["ec"].str.join("|")
    parity_sugars.rename(columns = {"score": "parityScore:float", "pdbl_subparity": "subParityScore:float", "ec": "ecID:string[]"}, inplace = True)

    parity_rels = pd.concat([parity_bl, parity_sugars])
    parity_rels.to_csv(f"{args.outdir}/bound_entity_parity_score_rels.csv.gz", compression = "gzip", sep = "\t", index = False)

    cath_bl_domain_ligand_interactions = cath_bl_domains[["cath_domain", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "bound_ligand_unique_id"]].drop_duplicates().rename(columns = {"bound_ligand_unique_id": ":END_ID(be-id)"})
    cath_sugar_domain_ligand_interactions = cath_sugar_domains[["cath_domain", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "sugar_id"]].drop_duplicates().rename(columns = {"sugar_id": ":END_ID(be-id)"})
    cath_domain_ligand_interactions = pd.concat([cath_bl_domain_ligand_interactions,cath_sugar_domain_ligand_interactions]).drop_duplicates()
    cath_domain_ligand_interactions.rename(columns = {"cath_domain": ":START_ID(cath-domain-ID)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc" , "domain_ownership" : "interactionMode"}, inplace = True)
    cath_domain_ligand_interactions.to_csv(f"{args.outdir}/cath_domain_ligand_interactions.csv", sep = "\t", index = False)

    scop_bl_domain_ligand_interactions = scop_bl_domains[["scop_id", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "bound_ligand_unique_id"]].drop_duplicates().rename(columns = {"bound_ligand_unique_id": ":END_ID(be-id)"})
    scop_sugar_domain_ligand_interactions = scop_sugar_domains[["scop_id", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "sugar_id"]].drop_duplicates().rename(columns = {"sugar_id": ":END_ID(be-id)"})
    scop_domain_ligand_interactions = pd.concat([scop_bl_domain_ligand_interactions,scop_sugar_domain_ligand_interactions]).drop_duplicates()
    scop_domain_ligand_interactions.rename(columns = {"scop_id": ":START_ID(scop-domain-id)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_ownership" : "interactionMode"}, inplace = True)
    scop_domain_ligand_interactions.to_csv(f"{args.outdir}/scop_domain_ligand_interactions.csv", sep = "\t", index = False)

    interpro_bl_domain_ligand_interactions = interpro_bl_domains[["interpro_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "bound_ligand_unique_id"]].drop_duplicates().rename(columns = {"bound_ligand_unique_id": ":END_ID(be-id)"})
    interpro_sugar_domain_ligand_interactions = interpro_sugar_domains[["interpro_accession", "domain_contact_counts", "domain_contact_perc", "domain_hbond_counts", "domain_hbond_perc", "domain_ownership", "sugar_id"]].drop_duplicates().rename(columns = {"sugar_id": ":END_ID(be-id)"})
    interpro_domain_ligand_interactions = pd.concat([interpro_bl_domain_ligand_interactions,interpro_sugar_domain_ligand_interactions]).drop_duplicates()
    interpro_domain_ligand_interactions.rename(columns = {"interpro_accession": ":START_ID(ipr-domain-ID)", "domain_contact_counts" : "domainContactCounts", "domain_contact_perc": "domainContactPerc", "domain_hbond_counts" : "domainHbondCounts", "domain_hbond_perc" : "domainHbondPerc", "domain_ownership" : "interactionMode"}, inplace = True)
    interpro_domain_ligand_interactions.to_csv(f"{args.outdir}/interpro_domain_ligand_interactions.csv", sep = "\t", index = False)

    scop_pdb_protein_rels = pd.concat([scop_bl_domains[["pdb_id", "protein_entity_id"]], scop_sugar_domains[["pdb_id", "protein_entity_id"]]]).drop_duplicates()
    cath_pdb_protein_rels = pd.concat([cath_bl_domains[["pdb_id", "protein_entity_id"]], cath_sugar_domains[["pdb_id", "protein_entity_id"]]]).drop_duplicates()
    interpro_pdb_protein_rels = pd.concat([interpro_bl_domains[["pdb_id", "protein_entity_id"]], interpro_sugar_domains[["pdb_id", "protein_entity_id"]]]).drop_duplicates()

    pdb_protein_rels = pd.concat([scop_pdb_protein_rels, cath_pdb_protein_rels, interpro_pdb_protein_rels]).drop_duplicates()
    pdb_protein_rels.rename(columns = {"protein_entity_id": ":START_ID(pdbp-id)", "pdb_id": ":END_ID(pdb-id)"}, inplace = True)
    pdb_protein_rels.to_csv(f"{args.outdir}/pdb_protein_rels.csv", sep = "\t", index = False)

    cath_protein_rels = pd.concat([cath_bl_domains[["cath_domain", "protein_entity_id"]], cath_sugar_domains[["cath_domain", "protein_entity_id"]]]).drop_duplicates().rename(columns = {"cath_domain": ":START_ID(cath-domain-ID)", "protein_entity_id":":END_ID(pdbp-id)"})
    scop_protein_rels = pd.concat([scop_bl_domains[["scop_id", "protein_entity_id"]], scop_sugar_domains[["scop_id", "protein_entity_id"]]]).drop_duplicates().rename(columns = {"scop_id": ":START_ID(scop-domain-id)", "protein_entity_id":":END_ID(pdbp-id)"})
    interpro_protein_rels = pd.concat([interpro_bl_domains[["interpro_accession", "protein_entity_id"]], interpro_sugar_domains[["interpro_accession", "protein_entity_id"]]]).drop_duplicates().rename(columns = {"interpro_accession": ":START_ID(ipr-domain-ID)", "protein_entity_id":":END_ID(pdbp-id)"})
    cath_protein_rels.to_csv(f"{args.outdir}/cath_protein_rels.csv", sep = "\t", index = False)
    scop_protein_rels.to_csv(f"{args.outdir}/scop_protein_rels.csv", sep = "\t", index = False)
    interpro_protein_rels.to_csv(f"{args.outdir}/interpro_protein_rels.csv", sep = "\t", index = False)