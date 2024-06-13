#!/usr/bin/env python

import pandas as pd
import re
import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--cath_domain_ownership', metavar='cath_domain_ownership', type=str,
                        help = "path to cath domain ownership file")
    parser.add_argument('--scop_domain_ownership', metavar='scop_domain_ownership', type=str,
                        help = "path to scop domain ownership file")
    parser.add_argument('--scop2_sf_domain_ownership', metavar='scop2_sf_domain_ownership', type=str,
                        help = "path to scop2_sf domain ownership file")
    parser.add_argument('--scop2_fa_domain_ownership', metavar='scop2_fa_domain_ownership', type=str,
                        help = "path to scop2_fa domain ownership file")
    parser.add_argument('--pfam_domain_ownership', metavar='pfam_domain_ownership', type=str,
                        help = "path to pfam domain ownership file")
    parser.add_argument('--superfamily_domain_ownership', metavar='superfamily_domain_ownership', type=str,
                        help = "path to superfamily domain ownership file")
    parser.add_argument('--gene3dsa_domain_ownership', metavar='gene3dsa_domain_ownership', type=str,
                        help = "path to gene3dsa domain ownership file")
    parser.add_argument('--scores_file', type=str, help='The scores file')
    parser.add_argument('--score_cutoff', type=float, default = 0.40, help='The score cutoff for cognate ligands')
    args = parser.parse_args()

    scores = pd.read_pickle(args.scores_file)
    scores_cognate = scores.loc[scores.score >= args.score_cutoff].copy()
    scores_cognate["ec"] = scores_cognate["ec"].str.split(",")
    scores_exploded = scores_cognate.explode("ec")
    scores_mask = scores_exploded.groupby(['ec', 'pdb_ligand'])['score'].transform(max) 
    scores_max = scores_exploded.loc[scores_exploded.score == scores_mask]
    
    cogligs = pd.read_pickle(args.cognate_ligands)
    cogligs.rename(columns = {"uniqueID": "cogliguid"}, inplace = True)
    
    core_cols = ['pdb_id', 'uniqueID', 'xref_db', 'xref_db_acc', 'xref_db_version', 'domain_accession', 'domain_type', 'derived_from', 'hetCode', 'type', 'bound_ligand_struct_asym_id','assembly_chain_id_ligand', 'assembly_chain_id_protein', 'proteinStructAsymID', 'bound_ligand_residue_interactions', 'domain_residue_interactions', 'domain_contact_counts', 'domain_hbond_counts', 'domain_covalent_counts', 'total_contact_counts', 'domain_contact_perc', 'num_non_minor_domains', 'domain_ownership', 'protein_entity_ec', 'ec_list', 'display_ec_list', 'chain_ec_list', 'pdb_ec_list', 'cognate_ligand', 'score', 'parity_match', 'parity_smarts', 'compound_name', 'ligand_db', 'compound_reaction', 'isCofactor']
    cath_cols = ['cath_code']
    scop_cols = ['sccs', 'domain_sunid']
    scop2_sf_cols = ['SCOPCLA']
    scop2_fa_cols = ['SCOPCLA']
    pfam_cols = ["pfam", "clan", "clan_description", "clan_comment", "pfam_description"]
    superfamily_cols = ["domain_description"]
    gene3dsa_cols = []

    #produce flat file where cognate ligands are present (bound entities without cognate ligands are identified from )
    cath_domains = pd.read_csv(f"{args.cath_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str", "cath_architecture": "str", "cath_class": "str", "cath_topology": "str", "cath_homologous_superfamily": "str"}, sep = "\t")
    cath_domains["pdb_ec_list"] = cath_domains["pdb_ec_list"].str.split(",")
    cath_exploded = cath_domains.explode("pdb_ec_list")
    cath_merged = cath_exploded.merge(scores_max, left_on = ["ligand_uniqueID", "ec_list"], right_on = ["pdb_ligand", "ec"], how = "inner")
    cath_merged_coglig = cath_merged.merge(cogligs, left_on = ["cognate_ligand","pdb_ec_list"], right_on = ["cogliguid", "entry"], how = "left")
    cath_merged_subset = cath_merged_coglig[core_cols + cath_cols]
    cath_merged_subset.to_csv("cath_cognate_ligands.csv.gz", compression = "gzip")

    scop_domains = pd.read_csv(f"{args.scop_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    scop_domains["pdb_ec_list"] = scop_domains["pdb_ec_list"].str.split(",")
    scop_exploded = scop_domains.explode("pdb_ec_list")
    scop_merged = scop_exploded.merge(scores_max, left_on = ["ligand_uniqueID", "ec_list"], right_on = ["pdb_ligand", "ec"], how = "inner")
    scop_merged_coglig = scop_merged.merge(cogligs, left_on = ["cognate_ligand","pdb_ec_list"], right_on = ["cogliguid", "entry"], how = "left")
    scop_merged_subset = scop_merged_coglig[core_cols + scop_cols]
    scop_merged_subset.to_csv("scop_cognate_ligands.csv.gz", compression = "gzip")
    
    pfam_domains = pd.read_csv(f"{args.pfam_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    superfamily_domains =pd.read_csv(f"{args.superfamily_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    gene3dsa_domains = pd.read_csv(f"{args.gene3dsa_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str", "cath_architecture": "str", "cath_class": "str", "cath_topology": "str", "cath_homologous_superfamily": "str"}, sep = "\t")
    scop2_sf_domains = pd.read_csv(f"{args.scop2_sf_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    scop2_fa_domains = pd.read_csv(f"{args.scop2_fa_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str"}, sep = "\t")
    

    #find all bound ligands where domains are all from the same chain, as contigs and af structures are not multichain (yet!) - data to be used for domain interaction tools.
    cath_assembly_chain_grouped = cath_domains.groupby('uniqueID')['chainUniqueID'].nunique().reset_index()
    cath_assembly_chain_grouped_filtered = cath_assembly_chain_grouped[cath_assembly_chain_grouped['chainUniqueID'] == 1]['uniqueID']
    cath_single_chain = cath_domains[cath_domains['uniqueID'].isin(cath_assembly_chain_grouped_filtered)].copy()
    cath_single_chain_exploded = cath_single_chain.explode("pdb_ec_list")
    cath_single_chain_merged = cath_single_chain_exploded.merge(scores_max, left_on = ["ligand_uniqueID", "ec_list"], right_on = ["pdb_ligand", "ec"], how = "left") #we use a left merge, and transfer pdb ligands. When there is a cognate match we are able to make this annotation, if not we can provide the pdb ligand. 
    cath_merged_nonminor = cath_single_chain_merged.loc[cath_single_chain_merged.domain_ownership != "minor"].copy()
    cath_merged_nonminor["combined_interaction"] = cath_merged_nonminor["cath_code"] + ":" + cath_merged_nonminor["domain_ownership"]
    cath_merged_nonminor["cath_segments_dict"] = cath_merged_nonminor.cath_segments_dict.apply(eval)
    cath_merged_nonminor["cath_min_start"] = cath_merged_nonminor["cath_segments_dict"].apply(lambda x: min([int(re.search(r"START=([-+\d]+)", y.get("SRANGE")).group(1)) for y in x]))
    cath_merged_nonminor_subset = cath_merged_nonminor[["pdb_id", "uniqueID", "hetCode", "pdb_ligand", "cognate_ligand", "combined_interaction"]].groupby(["pdb_id", "uniqueID", "hetCode", "pdb_ligand", "cognate_ligand"]).agg({"combined_interaction":list}).reset_index()#.head(10).T[10:]#[[""]]
    cath_merged_nonminor_subset["combined_interaction"] = cath_merged_nonminor_subset["combined_interaction"].str.join(";")
    cath_merged_nonminor_subset.to_csv("cath_single_chain_domain_interactions.tsv.gz", index = False, sep = "\t")

if __name__ == "__main__":
    main()
