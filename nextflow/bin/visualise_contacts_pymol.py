#!/usr/bin/env python

from pymol import cmd
import pandas as pd
import argparse
from pathlib import Path
import re
import math 

def visualise_contacts(contacts_df, domain_df):
    for i, contact in contacts_df.iterrows():
        obj_name = f"contact_{i+1}"
        atom1 = contact["atom1"]
        atom2 = contact["atom2"]
        if "," in atom1:
            atom1 = atom1.replace(",", "+")
            #create a psuedoatom at center of multiple atoms for line
            cmd.pseudoatom(f"{obj_name}_ps1", atom1)
            atom1 = f"{obj_name}_ps1"
        if "," in atom2:
            atom2 = atom2.replace(",", "+")
            #create a psuedoatom at center of multiple atoms for line
            cmd.pseudoatom(f"{obj_name}_ps2", atom2)
            atom2 = f"{obj_name}_ps2"
        #draw distance
        cmd.distance(obj_name, atom1, atom2)
        cmd.select(f"{obj_name}_sele", f"{atom1.rsplit('/', 1)[0]} + {atom2.rsplit('/', 1)[0]}")
        cmd.show("lines", f"{obj_name}_sele")
        cmd.hide("labels", obj_name)
        cmd.pseudoatom(f"{obj_name}_label", f"{atom1} + {atom2}", label=",".join(contact.contact))
        cmd.color("blue",f"{obj_name}") ##TODO - ALTER THIS COLOUR ACCORDING TO CONTACT TYPE ? 
    colour_list = ["blue", "red", "yellow", "pink", "black" , "orange"]
    colour_list = colour_list * math.ceil(domain_df.domain_accession.nunique()/len(colour_list))
    for i,domain in domain_df.reset_index(drop = True).iterrows():
        domain_name = domain["domain_accession"]
        cmd.select(f"{domain_name}", domain["selections"])
        cmd.color(colour_list[i], f"{domain_name}")
            

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--contacts_file', type=str, help='The contacts file')
    parser.add_argument('--structure_file', type=str, help='The structure file')
    parser.add_argument('--output_prefix', type=str, default = 'contacts', help='The output pse file')
    args = parser.parse_args()

    contacts = pd.read_csv(args.contacts_file)
    #struct name is assigned basename of pdb file
    struct_name = Path(re.sub(r".gz$", "", args.structure_file)).stem
    contacts["bgn"] = contacts["bgn"].apply(eval)
    contacts["end"] = contacts["end"].apply(eval)
    contacts["contact"] = contacts["contact"].apply(eval)
    contacts["auth_seq_range"] = contacts["auth_seq_range"].apply(eval)
    contacts["auth_seq_range"] = contacts["auth_seq_range"].apply(lambda x: "+".join([str(y).strip() for y in x])) #sometimes there is a space following last value in a multi-residue range

    contacts["atom1"] = contacts.apply(lambda x: f"/{struct_name}//{x['bgn'].get('auth_asym_id')}/{x['bgn'].get('auth_seq_id')}/{x['bgn'].get('auth_atom_id')}", axis = 1)
    contacts["atom2"] = contacts.apply(lambda x: f"/{struct_name}//{x['end'].get('auth_asym_id')}/{x['end'].get('auth_seq_id')}/{x['end'].get('auth_atom_id')}", axis = 1)

    for group, contact_group in contacts.groupby("xref_db"):
        domains = contact_group[["domain_accession", "auth_seq_range", "end_auth_asym_id"]].drop_duplicates(subset = ["domain_accession", "end_auth_asym_id"])
        domains["selections"] = f"/{struct_name}//" + domains["end_auth_asym_id"] + "/" + domains["auth_seq_range"]
        cmd.load(args.structure_file)
        visualise_contacts(contact_group, domains)
        cmd.save(f"{args.output_prefix}_{struct_name}_{group}.pse")
        cmd.delete("all")
    
    for ligand in contacts[["uniqueID", "bound_entity_pdb_residues"]].drop_duplicates():
        print(ligand)
        
if __name__ == "__main__":
    main()

#can use bound_entity_pdb_residues potentially to idnetify these in pymol