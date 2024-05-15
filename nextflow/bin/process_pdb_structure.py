#!/usr/bin/env python

from gemmi import cif
import pandas as pd
import argparse
from pathlib import Path
import sys
import re

def pattern_to_range(pattern):
    start, end = map(int, re.search(r'(\d+)-(\d+)', pattern).groups())
    return ",".join([str(x) for x in range(start, end + 1)])


def main():
    """
    This script processes the cif file of a structure and extracts the information needed to generate the arpeggio query file.
    The script will output a csv file containing the arpeggio query and a pickled dataframe containing information on the bound
    ligands in the structure.
    Script may be expected to fail if:
        1. The cif file does not have the right number of values in loop _atom_site (exitcode 99)
        2. The assembly contains more than 30 struct_asym ids (exitcode 100) - run these structures individually and combine manually.
        3. No domains are found in the cif file (exitcode 101)
        4. No bound entities are found in the cif file/assembly mapping (exitcode 102)

    """

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--cif', type=str, help='cif file containing PDBe updated structure, may be gzipped')
    parser.add_argument('--bio_h', type=str, help='cif file containing PDBe protonated structure, may be gzipped')
    parser.add_argument('--pdb_id', type=str, help='pdb id of the structure')
    parser.add_argument('--assembly_id', type=str, help='assembly id of the structure')
    args = parser.parse_args()

    try:
        doc = cif.read(args.cif)
    except ValueError as e:
        # know of at least one case of this error: 8j07_updated, so we catch it here and allow it in nextflow pipelien by specific exit code
        if "Wrong number of values in loop _atom_site" in str(e):
            sys.exit(100)
        else:
            sys.exit(1)

    block = doc.sole_block()

    #we read in the protonated structure to verify the mapping between updated and protonated structure
    #e.g. _ or - separator in oper id (see 2c16 for example)

    bioh_doc = cif.read(args.bio_h)
    bioh_block = bioh_doc.sole_block()

    bioh_struct_asym = pd.DataFrame(bioh_block.find(['_struct_asym.id']), columns = ["id"])
    bioh_struct_asym["sep"] = bioh_struct_asym.id.str.extract("^.+([-|_]).+$")
    separators = bioh_struct_asym["sep"].dropna().unique()
    if len(separators) == 0:
        print("No separators detected in assembly struct_asym")
        separator = ""
    elif len(separators) == 1:
        separator = separators[0]
    else:
        print("Mixed separators detected, exiting")
        sys.exit(1)

    entity_info = pd.DataFrame(block.find(['_entity.id', '_entity.pdbx_description', '_entity.formula_weight']), columns = ["entity_id", "description", "molweight"])
    entity_info["description"] = entity_info["description"].str.strip("\"|'")
    ##see this webinar for details https://pdbeurope.github.io/api-webinars/webinars/web5/arpeggio.html
    assembly_info = pd.DataFrame(block.find(['_pdbx_struct_assembly_gen.assembly_id', '_pdbx_struct_assembly_gen.oper_expression', '_pdbx_struct_assembly_gen.asym_id_list']), columns = ["assembly_id", "oper_expression", "asym_id_list"])
    assembly_info = assembly_info.loc[assembly_info.assembly_id == args.assembly_id].copy() #preferred assembly id
    #some structures have a range in the format '(1-60)' - expand this before splitting, see 1al0 for example, some structures have a range in the format 1-60, see 6rjf for example
    #so, first strip all brackets from the oper expression and then expand the range
    assembly_info["oper_expression"] = assembly_info["oper_expression"].str.strip("()'")
    assembly_info.loc[assembly_info["oper_expression"].str.match("\d+-\d+"), "oper_expression"] = assembly_info.loc[assembly_info["oper_expression"].str.match("\d+-\d+"), "oper_expression"].apply(pattern_to_range)
    #observe some ; and \n in asym_id_list (see 3hye for example) -  so strip ; and \n; from start and end of string before splitting - will expand this if necessary on more errors
    assembly_info["oper_expression"] = assembly_info["oper_expression"].str.strip("\n;")
    assembly_info["oper_expression"] = assembly_info["oper_expression"].str.split(",")
    assembly_info["asym_id_list"] = assembly_info["asym_id_list"].str.strip("\n;").str.split(",") #asym id here is the struct
    assembly_info_exploded = assembly_info.explode("oper_expression").explode("asym_id_list").rename(columns = {"asym_id_list": "struct_asym_id"})
    #the oper_expression for the identity operation does not receive an _[digit] suffix, so we need to remove the oper_expression from it - see 2r9p for example or 4cr7
    #this is necessary for matching the pdb-h format
    oper_list = pd.DataFrame(block.find(['_pdbx_struct_oper_list.id', '_pdbx_struct_oper_list.type']), columns = ["oper_expression_id", "type"])

    assembly_info_exploded_oper = assembly_info_exploded.merge(oper_list, left_on = "oper_expression", right_on = "oper_expression_id", how = "left", indicator = True)
    assert(len(assembly_info_exploded_oper.loc[assembly_info_exploded_oper._merge != "both"]) == 0)
    assembly_info_exploded_oper.loc[assembly_info_exploded_oper.type == "'identity operation'", "oper_expression"] = ""
    assembly_info_exploded_oper.drop(columns = ["oper_expression_id", "type" , "_merge"], inplace = True)

    #arpeggio can spend over 30 hours running on some structures which have a high molecular weight like capsid structures. In order to build the database in a realistic time, we filter out structures with assemblies > Xkda
    struct_asym_info = pd.DataFrame(block.find(['_struct_asym.id', '_struct_asym.entity_id']), columns = ["struct_asym_id", "entity_id"])
    struct_asym_info = struct_asym_info.merge(entity_info, on = "entity_id", indicator = True)
    assert(len(struct_asym_info.loc[struct_asym_info._merge != "both"]) == 0)
    struct_asym_info.drop(columns = "_merge", inplace = True)

    struct_asym_info = struct_asym_info.merge(assembly_info_exploded_oper, on = "struct_asym_id", indicator = True)
    assert(len(struct_asym_info.loc[struct_asym_info._merge != "both"]) == 0)
    #unknown atom or ion can result in ? in molweight - e.g. 1svu, so filter the df before taking molweight sum (also filter dot as we know these can occur elsewhere)
    assembly_molwt_kda = struct_asym_info.loc[struct_asym_info.molweight.isin(["?", "."]) == False].molweight.astype("float").sum() / 1000

    if assembly_molwt_kda >= 250:
        print("Large structure detected, run individually instead of in pipeline. exiting")
        sys.exit(100)

    ##check if the structure has any domains before continuing (if not the structure should not be processed) - e.g. 1lg1 (1i8q shows why auth doesn't work)
    domain_info_dataframe = pd.DataFrame(block.find("_pdbx_sifts_xref_db_segments.", ["entity_id", "asym_id", "xref_db", "xref_db_acc", "domain_name", "segment_id", "instance_id", "seq_id_start", "seq_id_end"]), 
        columns = ["entity_id", "asym_id", "xref_db", "xref_db_acc", "domain_name", "segment_id", "instance_id", "seq_id_start", "seq_id_end"])
    ##for now, we are ready to handle annotations from CATH, SCOP, SCOP2B, Pfam and InterPro, so we filter to this
    domain_info_dataframe_filtered = domain_info_dataframe.loc[domain_info_dataframe.xref_db.isin(["CATH", "SCOP", "SCOP2B", "Pfam", "InterPro"])]
    if len(domain_info_dataframe_filtered) == 0:
        print("No domains found in cif file [CATH, SCOP, SCOP2B, Pfam, InterPro]")
        sys.exit(101)
    #switched auth_asym_id to pdb_asym_id (which is a pointer to atom_site auth asym id and seems to hold up better for branch structures in the mapping to pdb-h structure)
    branched_seq_info = pd.DataFrame(block.find(['_pdbx_branch_scheme.asym_id', '_pdbx_branch_scheme.mon_id', '_pdbx_branch_scheme.entity_id', '_pdbx_branch_scheme.pdb_seq_num', '_pdbx_branch_scheme.pdb_asym_id', '_pdbx_branch_scheme.auth_seq_num']), columns = ["bound_ligand_struct_asym_id", "hetCode", "entity_id", "pdb_seq_num", "auth_asym_id", "auth_seq_num"])
    branched_seq_info_merged =  pd.DataFrame([], columns = ['bound_ligand_struct_asym_id', 'hetCode', 'entity_id', 'pdb_seq_num', 'auth_asym_id', 'auth_seq_num', 'descriptor'])
    if len(branched_seq_info) > 0:
        branched_seq_info["hetCode"] = "SUGAR"
        branched_sugar_info = pd.DataFrame(block.find(['_pdbx_entity_branch_descriptor.entity_id', '_pdbx_entity_branch_descriptor.descriptor', '_pdbx_entity_branch_descriptor.type']), columns = ["entity_id", "descriptor", "type"])
        branched_sugar_info_wurcs = branched_sugar_info.loc[branched_sugar_info.type == "WURCS"].groupby("entity_id").first()
        if len(branched_sugar_info_wurcs) == 0:
            print("No sugars with WURCS representation found")
        branched_sugar_info_wurcs["descriptor"] = branched_sugar_info_wurcs["descriptor"].str.strip("\"|'")
        branched_seq_info_merged = branched_seq_info.merge(branched_sugar_info_wurcs, on = "entity_id", indicator = True)
        assert(len(branched_seq_info_merged.loc[branched_seq_info_merged._merge != "both"]) == 0)
        branched_seq_info_merged.drop(columns = ["_merge", "type"], inplace = True)
        branched_seq_info_merged["type"] = "sugar"
        branched_seq_info_merged["pdb_ins_code"] = "" #add blank ins code for branch structures for downstream agg

    nonpoly_info = pd.DataFrame(block.find(['_pdbx_entity_nonpoly.entity_id', '_pdbx_entity_nonpoly.comp_id']), columns = ["entity_id", "hetCode"])
    nonpoly_info = nonpoly_info.loc[nonpoly_info.hetCode.isin(["HOH", "UNL", "UNX"]) == False]
    nonpoly_seq_info_filtered = pd.DataFrame([], columns = ['bound_ligand_struct_asym_id', 'entity_id', 'pdb_seq_num', 'auth_asym_id', 'auth_seq_num', 'hetCode'])
    if len(nonpoly_info) > 0:
        nonpoly_seq_info = pd.DataFrame(block.find(['_pdbx_nonpoly_scheme.asym_id', '_pdbx_nonpoly_scheme.entity_id', '_pdbx_nonpoly_scheme.pdb_seq_num', '_pdbx_nonpoly_scheme.pdb_strand_id', '_pdbx_nonpoly_scheme.auth_seq_num', '_pdbx_nonpoly_scheme.pdb_ins_code']), columns = ["bound_ligand_struct_asym_id","entity_id", "pdb_seq_num", "auth_asym_id", "auth_seq_num", "pdb_ins_code"])
        nonpoly_seq_info_merged = nonpoly_seq_info.merge(nonpoly_info, on = "entity_id", indicator = True)
        assert(len(nonpoly_seq_info_merged.loc[nonpoly_seq_info_merged._merge != "both"])  == 0)
        nonpoly_seq_info_merged.drop(columns = "_merge", inplace = True)
        nonpoly_seq_info_merged["pdb_ins_code"] = nonpoly_seq_info_merged["pdb_ins_code"].fillna("").str.replace("\?|\.", "",regex = True)
        nonpoly_seq_info_filtered = nonpoly_seq_info_merged.loc[nonpoly_seq_info_merged.hetCode.isin(["HOH", "UNL"]) == False]
        nonpoly_seq_info_filtered["type"] = "ligand"

    if len(branched_seq_info_merged) > 0 or len(nonpoly_seq_info_filtered) > 0:
        bound_entity_info = pd.concat([branched_seq_info_merged, nonpoly_seq_info_filtered])
        bound_entity_info = bound_entity_info.merge(entity_info, on = "entity_id", how = "left")
        bound_entity_info_assembly = bound_entity_info.merge(assembly_info_exploded_oper, left_on = "bound_ligand_struct_asym_id", right_on = "struct_asym_id", how = "inner") #inner join only keep entities in the assembly
        #assembly may not include all bound entities. e.g. 3m43 GOL not in assembly
        if len(bound_entity_info_assembly) == 0:
            print("No bound entities mapped to the assembly operations for preferred assembly")
            sys.exit(102)
        #assert(len(bound_entity_info_assembly.loc[bound_entity_info_assembly._merge != "both"]) == 0)
        #originally beleived the protonated structure for branch used struct asym id, but it seems to use auth asym id so following two lines are identical and could be deduplicated once we are sure this is true
        bound_entity_info_assembly.loc[(bound_entity_info_assembly.type == "ligand"), "assembly_chain_id_ligand"] = bound_entity_info_assembly.loc[(bound_entity_info_assembly.type == "ligand"), "auth_asym_id"] + separator + bound_entity_info_assembly.loc[(bound_entity_info_assembly.type == "ligand"), "oper_expression"]
        bound_entity_info_assembly.loc[(bound_entity_info_assembly.type == "sugar"), "assembly_chain_id_ligand"] = bound_entity_info_assembly.loc[(bound_entity_info_assembly.type == "sugar"), "auth_asym_id"] + separator + bound_entity_info_assembly.loc[(bound_entity_info_assembly.type == "sugar"), "oper_expression"]
        bound_entity_info_assembly["assembly_chain_id_ligand"] = bound_entity_info_assembly["assembly_chain_id_ligand"].str.strip(separator)
        bound_entity_info_assembly.drop(columns = ["struct_asym_id", "oper_expression"], inplace = True) #"_merge", 

        bound_entity_info_grouped = bound_entity_info_assembly.groupby(["bound_ligand_struct_asym_id", "assembly_chain_id_ligand", "entity_id"]).agg({"pdb_ins_code": "first", "hetCode": "first", "descriptor": "first", "description": "first", "type": "first", "auth_seq_num": list, "pdb_seq_num": list}).reset_index()
        bound_entity_info_grouped["bound_molecule_display_id"] = bound_entity_info_grouped.groupby(["assembly_chain_id_ligand", "entity_id"]).ngroup().transform(lambda x: "bm"+ str(x+1)) #we could sort here to try and put bm of identical entities together ? or actually we may want to not sort the groupby as it already is
        bound_entity_info_grouped["uniqueID"] = args.pdb_id + "_" + bound_entity_info_grouped["bound_molecule_display_id"] + "_" + bound_entity_info_grouped["bound_ligand_struct_asym_id"]
        bound_entity_info_grouped.loc[bound_entity_info_grouped.type == "ligand", "arpeggio"] = bound_entity_info_grouped.loc[bound_entity_info_grouped.type == "ligand"].apply(lambda x: [f"/{x['assembly_chain_id_ligand']}/{y}{x['pdb_ins_code']}/" for y in x["pdb_seq_num"]], axis = 1)
        bound_entity_info_grouped.loc[bound_entity_info_grouped.type == "sugar", "arpeggio"] = bound_entity_info_grouped.loc[bound_entity_info_grouped.type == "sugar"].apply(lambda x: [f"/{x['assembly_chain_id_ligand']}/{y}/" for y in x["pdb_seq_num"]], axis = 1)
        bound_entity_info_grouped["bound_entity_auth_residues"] = bound_entity_info_grouped["auth_seq_num"].str.join("|")
        bound_entity_info_grouped["pdb_seq_num"] = bound_entity_info_grouped["pdb_seq_num"].str.join("|")
        bound_entity_info_grouped["entity_id"] = bound_entity_info_grouped["entity_id"].astype(int)
        bound_entity_info_grouped.rename(columns = {"pdb_seq_num" : "bound_entity_pdb_residues", "entity_id": "ligand_entity_id_numerical"}, inplace = True)
        bound_entity_info_grouped["arpeggio"].explode().to_csv(f"{args.pdb_id}_arpeggio.csv", index = False, header = None)
        bound_entity_info_grouped.to_pickle(f"{args.pdb_id}_bound_entity_info.pkl")
    else:
        print("No bound entities found in cif file")
        sys.exit(102)

if __name__ == "__main__":
    main()
