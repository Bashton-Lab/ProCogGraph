#!/usr/bin/env python
from pdbecif.mmcif_io import CifFileReader
from Bio.ExPASy import Enzyme as EEnzyme
from pathlib import Path
import requests
import os
from pprint import pprint # for later pretty printing only

#a limitation of this function is that it selects only one representative EC number for transferred entries
def get_terminal_record(entry, row, df):
    entry = row.ID
    pattern = r"[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+"
    while row.DE.startswith("Transferred entry: "):
        transfers = re.findall(pattern, row.DE)
        #when updating, if multiple possible transfers, selecting only first.
        row = df.loc[df.ID == transfers[0]].iloc[0]
    return row.ID

import json
from urllib.parse import quote
def get_glycoct_from_wurcs(wurcs):
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

def get_csdb_from_glycoct(glycoct):
    if glycoct is np.nan:
        return np.nan
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


import numpy as np
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
    ec_list = ec_records_df.ID.unique()
    
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
    return(df_merged)

from neo4j import __version__ as neo4j_version,  GraphDatabase
print(f"Neo4j python package version: {neo4j_version}")
#class is from https://towardsdatascience.com/neo4j-cypher-python-7a919a372be7
import pandas as pd
class Neo4jConnection:
    
    def __init__(self, uri, user, pwd):
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            print("Failed to create the driver:", e)
        
    def close(self):
        if self.__driver is not None:
            self.__driver.close()
        
    def query(self, query, db=None):
        assert self.__driver is not None, "Driver not initialized!"
        session = None
        response = None
        try: 
            session = self.__driver.session(database=db) if db is not None else self.__driver.session() 
            response = list(session.run(query))
        except Exception as e:
            print("Query failed:", e)
        finally: 
            if session is not None:
                session.close()
        return response

def parse_table_data(elem):
    data = {}
    for row in elem.findall('row'):
        row_data = {}
        for field in row.findall('field'):
            row_data[field.attrib['name']] = field.text
        data[len(data)] = row_data
    return data

import pandas as pd
import xml.etree.ElementTree as ET
import re
import argparse

from bs4 import BeautifulSoup
from urllib.parse import quote
def get_smiles_from_csdb(csdb_linear):
    if csdb_linear is np.nan:
        return np.nan
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

def main():

    parser = argparse.ArgumentParser(description = 'TO DO')
    parser.add_argument('--neo4j_bolt_uri', default = 'bolt://localhost:7687', type = str,
        help = "")
    parser.add_argument('--neo4j_user', default = 'neo4j', type = str,
        help = "")
    parser.add_argument('--neo4j_password', type = str,
        help = "")
    parser.add_argument('--outdir', type = str, default = "parity_calcs",
        help = ""),
    parser.add_argument('--enzyme_dat_file', type = str, default = "enzyme.dat",
        help = ""),
    
    args = parser.parse_args()

    print("Connecting to neo4j")
    conn = Neo4jConnection(uri=f"{args.neo4j_bolt_uri}", user=f"{args.neo4j_user}", pwd=f"{args.neo4j_password}")
    print("Connected to neo4j")
    print("Generating EC record dataframe")
    with open(f"{args.enzyme_dat_file}") as handle:
        ec_records = EEnzyme.parse(handle)
        ec_records_list = []
        for record in ec_records: 
            ec_record_series = pd.Series(record)
            ec_records_list.append(ec_record_series)


    ec_records_df = pd.DataFrame(ec_records_list)
    ec_records_df["TRANSFER"] = ec_records_df.apply(lambda x: get_terminal_record(x["ID"], x, ec_records_df), axis = 1)
    ec_records_df["TRANSFER"] = ec_records_df["TRANSFER"].fillna(ec_records_df.ID)

    #how do we store the neo4j cypher queries? in a separate file? yaml?
    #should consider how we can use the symmetry operator information to condense this down? multiple bound molecules can conist of the same bound ligands, differentiated by SYM_OPERATOR
    #in the bl instances, we are not making use of the symmetry operator information. Or is it covered by having the bm info there?
    
    cath_pdb_residue_interactions_query_distinct_bl = """
    MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity  {TYPE : 'b', POLYMER_TYPE : 'B'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[cd:IS_IN_CATH_DOMAIN]->(c:CATH),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity  {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']" 
    AND cd.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    p.UNIQID as ligand_entity_id, 
    p.ID as ligand_entity_id_numerical,
    p.DESCRIPTION as ligand_entity_description,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    cd.AUTH_ASYM_ID as cath_chain_id,
    c.DOMAIN as cath_domain,
    c.CLASS as cath_class,
    c.ARCH as cath_architecture,
    c.TOPOL as cath_topology,
    c.HOMOL as cath_homology,
    c.NAME as cath_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.CHEM_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    cath_pdb_residue_interactions_query_distinct_sugar = """
    MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity  {TYPE : 's', POLYMER_TYPE : 'S'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[cd:IS_IN_CATH_DOMAIN]->(c:CATH),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity  {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']" 
    AND cd.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    p.UNIQID as ligand_entity_id, 
    p.ID as ligand_entity_id_numerical,
    p.DESCRIPTION as ligand_entity_description,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    cd.AUTH_ASYM_ID as cath_chain_id,
    c.DOMAIN as cath_domain,
    c.CLASS as cath_class,
    c.ARCH as cath_architecture,
    c.TOPOL as cath_topology,
    c.HOMOL as cath_homology,
    c.NAME as cath_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.CHEM_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    scop_pdb_residue_interactions_query_distinct_bl = """
    MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 'b', POLYMER_TYPE : 'B'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[sd:IS_IN_SCOP_DOMAIN]->(s:SCOP),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']" 
    AND sd.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id,
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    s.SUNID as scop_sunid,
    s.DESCRIPTION as scop_description,
    s.SCCS as scop_sccs,
    sd.CLASS_ID as scop_class_id,
    sd.FOLD_ID as scop_fold_id,
    sd.SUPERFAMILY_ID as scop_superfamily_id,
    sd.SCOP_ID as scop_id,
    sd.AUTH_ASYM_ID as scop_chain_id,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    scop_pdb_residue_interactions_query_distinct_sugar = """
    MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 's', POLYMER_TYPE : 'S'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[sd:IS_IN_SCOP_DOMAIN]->(s:SCOP),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry) 
    WHERE
    e.EC IS NOT NULL AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']" 
    AND sd.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id, 
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    s.SUNID as scop_sunid,
    s.DESCRIPTION as scop_description,
    s.SCCS as scop_sccs,
    sd.CLASS_ID as scop_class_id,
    sd.FOLD_ID as scop_fold_id,
    sd.SUPERFAMILY_ID as scop_superfamily_id,
    sd.SCOP_ID as scop_id,
    sd.AUTH_ASYM_ID as scop_chain_id,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2"""

    interpro_pdb_residue_interactions_query_distinct_bl_d = """MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 'b', POLYMER_TYPE : 'B'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[ip:IS_IN_INTERPRO]->(i:Interpro),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND i.ENTRY_TYPE = "D" AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']"
    AND ip.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id, 
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    i.ABBREV as interpro_type,
    i.INTERPRO_ACCESSION as interpro_accession,
    i.NAME as interpro_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    interpro_pdb_residue_interactions_query_distinct_sugar_d = """MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 's', POLYMER_TYPE : 'S'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[ip:IS_IN_INTERPRO]->(i:Interpro),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND i.ENTRY_TYPE = "D" AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']"
    AND ip.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id, 
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    i.ABBREV as interpro_type,
    i.INTERPRO_ACCESSION as interpro_accession,
    i.NAME as interpro_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    interpro_pdb_residue_interactions_query_distinct_bl_f = """MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 'b', POLYMER_TYPE : 'B'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[ip:IS_IN_INTERPRO]->(i:Interpro),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND i.ENTRY_TYPE = "F" AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']"
    AND ip.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id, 
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    i.ABBREV as interpro_type,
    i.INTERPRO_ACCESSION as interpro_accession,
    i.NAME as interpro_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    interpro_pdb_residue_interactions_query_distinct_sugar_f = """MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 's', POLYMER_TYPE : 'S'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[ip:IS_IN_INTERPRO]->(i:Interpro),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND i.ENTRY_TYPE = "F" AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']"
    AND ip.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id, 
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    i.ABBREV as interpro_type,
    i.INTERPRO_ACCESSION as interpro_accession,
    i.NAME as interpro_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    interpro_pdb_residue_interactions_query_distinct_bl_h = """MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 'b', POLYMER_TYPE : 'B'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[ip:IS_IN_INTERPRO]->(i:Interpro),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND i.ENTRY_TYPE = "H" AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']"
    AND ip.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id, 
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    i.ABBREV as interpro_type,
    i.INTERPRO_ACCESSION as interpro_accession,
    i.NAME as interpro_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """

    interpro_pdb_residue_interactions_query_distinct_sugar_h = """MATCH
    (bl:BoundLigand)<-[x:IS_AN_INSTANCE_OF]->(p:Entity {TYPE : 's', POLYMER_TYPE : 'S'})<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)-[Arp:HAS_ARP_CONTACT]->(pr:PDBResidue)-[ip:IS_IN_INTERPRO]->(i:Interpro),
    (pr)<-[y:HAS_PDB_RESIDUE]-(e:Entity {TYPE : 'p', POLYMER_TYPE : 'P'})<-[z:HAS_ENTITY]-(a:Entry)
    WHERE
    e.EC IS NOT NULL AND i.ENTRY_TYPE = "H" AND
    Arp.CONTACT_TYPE <> "['vdw_clash']" AND Arp.CONTACT_TYPE <> "['clash']"
    AND ip.STRUCT_ASYM_ID = Arp.STRUCT_ASYM_ID_2 AND
    NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND Arp.SYM_OP_1 = part.SYM_OPERATOR 
    AND bl.CHEM_COMP_ID in split(p.CHEM_COMP_LIST, ",")
    RETURN DISTINCT
    a.ID as pdb_id, 
    e.UNIQID as protein_entity_id,
    e.BEST_CHAIN_ID as pdb_best_chain_id,
    e.EC as protein_entity_ec,
    bl.UNIQID as bound_ligand_entity_id,
    p.UNIQID as ligand_entity_id, 
    p.DESCRIPTION as ligand_entity_description,
    p.ID as ligand_entity_id_numerical,
    p.TYPE as ligand_entity_type,
    p.BEST_CHAIN_ID as bound_ligand_best_chain_id, 
    p.POLYMER_TYPE as ligand_entity_polymer_type, 
    pr.UNIQID as pdb_residue_id,
    pr.CHEM_COMP_ID as pdb_residue_type,
    i.ABBREV as interpro_type,
    i.INTERPRO_ACCESSION as interpro_accession,
    i.NAME as interpro_name,
    bm.UNIQID as bound_molecule_id,
    bl.UNIQID as bound_ligand_id, 
    bl.AUTH_COMP_ID as bound_ligand_name,
    Arp.AUTH_SEQ_ID_2 as pdb_residue_auth_id,
    Arp.AUTH_SEQ_ID_1 as bound_ligand_auth_id,
    Arp.CONTACT_TYPE as contact_type, 
    Arp.DISTANCE as contact_distance, 
    Arp.INTERACTION_TYPE as interaction_type, 
    Arp.ATOM_1 as atom_1, 
    Arp.ATOM_2 as atom_2
    """
    Path(f"{args.outdir}").mkdir(parents=True, exist_ok=True)
    if not os.path.exists(f"{args.outdir}/cath_pdb_residue_interactions_distinct_bl.csv.gz"):
        print("Retrieving cath_pdb_residue_interactions_distinct_bl")
        cath_pdb_residue_interactions_distinct_bl = pd.DataFrame([dict(_) for _ in conn.query(cath_pdb_residue_interactions_query_distinct_bl, db='neo4j')])
        cath_pdb_residue_interactions_distinct_bl.to_csv(f"{args.outdir}/cath_pdb_residue_interactions_distinct_bl.csv.gz", compression = "gzip")
    else:
        print("Loading cath_pdb_residue_interactions_distinct_bl")
        cath_pdb_residue_interactions_distinct_bl = pd.read_csv(f"{args.outdir}/cath_pdb_residue_interactions_distinct_bl.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

    if not os.path.exists(f"{args.outdir}/cath_pdb_residue_interactions_distinct_sugar.csv.gz"):
        print("Retrieving cath_pdb_residue_interactions_distinct_sugar")
        cath_pdb_residue_interactions_distinct_sugar = pd.DataFrame([dict(_) for _ in conn.query(cath_pdb_residue_interactions_query_distinct_sugar, db='neo4j')])
        cath_pdb_residue_interactions_distinct_sugar.to_csv(f"{args.outdir}/cath_pdb_residue_interactions_distinct_sugar.csv.gz", compression = "gzip")
    else:
        print("Loading cath_pdb_residue_interactions_distinct_sugar")
        cath_pdb_residue_interactions_distinct_sugar = pd.read_csv(f"{args.outdir}/cath_pdb_residue_interactions_distinct_sugar.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

    if not os.path.exists(f"{args.outdir}/scop_pdb_residue_interactions_distinct_bl.csv.gz"):
        print("Retrieving scop_pdb_residue_interactions_distinct_bl")
        scop_pdb_residue_interactions_distinct_bl = pd.DataFrame([dict(_) for _ in conn.query(scop_pdb_residue_interactions_query_distinct_bl, db='neo4j')])
        scop_pdb_residue_interactions_distinct_bl.to_csv(f"{args.outdir}/scop_pdb_residue_interactions_distinct_bl.csv.gz", compression = "gzip")
    else:
        print("Loading scop_pdb_residue_interactions_distinct_bl")
        scop_pdb_residue_interactions_distinct_bl = pd.read_csv(f"{args.outdir}/scop_pdb_residue_interactions_distinct_bl.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    
    if not os.path.exists(f"{args.outdir}/scop_pdb_residue_interactions_distinct_sugar.csv.gz"):
        print("Retrieving scop_pdb_residue_interactions_distinct_sugar")
        scop_pdb_residue_interactions_distinct_sugar = pd.DataFrame([dict(_) for _ in conn.query(scop_pdb_residue_interactions_query_distinct_sugar, db='neo4j')])
        scop_pdb_residue_interactions_distinct_sugar.to_csv(f"{args.outdir}/scop_pdb_residue_interactions_distinct_sugar.csv.gz", compression = "gzip")
    else:
        print("Loading scop_pdb_residue_interactions_distinct_sugar")
        scop_pdb_residue_interactions_distinct_sugar = pd.read_csv(f"{args.outdir}/scop_pdb_residue_interactions_distinct_sugar.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    
    if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl.csv.gz"):
        if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_d.csv.gz"):
            print("Retrieving interpro_pdb_residue_interactions_distinct_bl_d")
            interpro_pdb_residue_interactions_distinct_bl_d = pd.DataFrame([dict(_) for _ in conn.query(interpro_pdb_residue_interactions_query_distinct_bl_d, db='neo4j')])
            interpro_pdb_residue_interactions_distinct_bl_d.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_d.csv.gz", compression = "gzip")
        else:
            print("Loading interpro_pdb_residue_interactions_distinct_bl_d")
            interpro_pdb_residue_interactions_distinct_bl_d = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_d.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        
        if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_f.csv.gz"):
            print("Retrieving interpro_pdb_residue_interactions_distinct_bl_f")
            interpro_pdb_residue_interactions_distinct_bl_f = pd.DataFrame([dict(_) for _ in conn.query(interpro_pdb_residue_interactions_query_distinct_bl_f, db='neo4j')])
            interpro_pdb_residue_interactions_distinct_bl_f.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_f.csv.gz", compression = "gzip")
        else:
            print("Loading interpro_pdb_residue_interactions_distinct_bl_f")
            interpro_pdb_residue_interactions_distinct_bl_f = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_f.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

        if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_h.csv.gz"):
            print("Retrieving interpro_pdb_residue_interactions_distinct_bl_h")
            interpro_pdb_residue_interactions_distinct_bl_h = pd.DataFrame([dict(_) for _ in conn.query(interpro_pdb_residue_interactions_query_distinct_bl_h, db='neo4j')])
            interpro_pdb_residue_interactions_distinct_bl_h.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_h.csv.gz", compression = "gzip")
        else:
            print("Loading interpro_pdb_residue_interactions_distinct_bl_h")
            interpro_pdb_residue_interactions_distinct_bl_h = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_h.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        interpro_pdb_residue_interactions_distinct_bl = pd.concat([interpro_pdb_residue_interactions_distinct_bl_d,
                                                                interpro_pdb_residue_interactions_distinct_bl_f,
                                                                interpro_pdb_residue_interactions_distinct_bl_h])
        interpro_pdb_residue_interactions_distinct_bl.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl.csv.gz", compression = "gzip")
    else:
        interpro_pdb_residue_interactions_distinct_bl = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    
    if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar.csv.gz"):
        if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_d.csv.gz"):
            print("Retrieving interpro_pdb_residue_interactions_distinct_sugar_d")
            interpro_pdb_residue_interactions_distinct_sugar_d = pd.DataFrame([dict(_) for _ in conn.query(interpro_pdb_residue_interactions_query_distinct_sugar_d, db='neo4j')])
            interpro_pdb_residue_interactions_distinct_sugar_d.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_d.csv.gz", compression = "gzip")
        else:
            print("Loading interpro_pdb_residue_interactions_distinct_sugar_d")
            interpro_pdb_residue_interactions_distinct_sugar_d = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_d.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

        if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_f.csv.gz"):
            print("Retrieving interpro_pdb_residue_interactions_distinct_sugar_f")
            interpro_pdb_residue_interactions_distinct_sugar_f = pd.DataFrame([dict(_) for _ in conn.query(interpro_pdb_residue_interactions_query_distinct_sugar_f, db='neo4j')])
            interpro_pdb_residue_interactions_distinct_sugar_f.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_f.csv.gz", compression = "gzip")
        else:
            print("Loading interpro_pdb_residue_interactions_distinct_sugar_f")
            interpro_pdb_residue_interactions_distinct_sugar_f = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_f.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
        
        if not os.path.exists(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_h.csv.gz"):
            print("Retrieving interpro_pdb_residue_interactions_distinct_sugar_h")
            interpro_pdb_residue_interactions_distinct_sugar_h = pd.DataFrame([dict(_) for _ in conn.query(interpro_pdb_residue_interactions_query_distinct_sugar_h, db='neo4j')])
            interpro_pdb_residue_interactions_distinct_sugar_h.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_h.csv.gz", compression = "gzip")
        else:
            print("Loading interpro_pdb_residue_interactions_distinct_sugar_h")
            interpro_pdb_residue_interactions_distinct_sugar_h = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_h.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)
    
        interpro_pdb_residue_interactions_distinct_sugar = pd.concat([interpro_pdb_residue_interactions_distinct_sugar_d,
                                                                interpro_pdb_residue_interactions_distinct_sugar_f,
                                                                interpro_pdb_residue_interactions_distinct_sugar_h])
        interpro_pdb_residue_interactions_distinct_sugar.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar.csv.gz", compression = "gzip")
    else:
        print("Loading interpro_pdb_residue_interactions_distinct_sugar")
        interpro_pdb_residue_interactions_distinct_bl = pd.read_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

    if not os.path.exists(f"{args.outdir}/bound_molecules_ligands"):
        print("Updating ligand EC records")
        cath_pdb_residue_interactions_distinct_bl_ec = get_updated_enzyme_records(cath_pdb_residue_interactions_distinct_bl, ec_records_df)
        scop_pdb_residue_interactions_distinct_bl_ec = get_updated_enzyme_records(scop_pdb_residue_interactions_distinct_bl, ec_records_df)
        interpro_pdb_residue_interactions_distinct_bl_ec = get_updated_enzyme_records(interpro_pdb_residue_interactions_distinct_bl, ec_records_df)

        cath_pdb_residue_interactions_distinct_bl_ec["uniqueID"] = cath_pdb_residue_interactions_distinct_bl_ec["bound_molecule_id"] + "_" + cath_pdb_residue_interactions_distinct_bl_ec["bound_ligand_id"]
        cath_pdb_residue_interactions_distinct_bl_ec["type"] = "ligand"
        scop_pdb_residue_interactions_distinct_bl_ec["uniqueID"] = scop_pdb_residue_interactions_distinct_bl_ec["bound_molecule_id"] + "_" + scop_pdb_residue_interactions_distinct_bl_ec["bound_ligand_id"]
        scop_pdb_residue_interactions_distinct_bl_ec["type"] = "ligand"
        interpro_pdb_residue_interactions_distinct_bl_ec["uniqueID"] = interpro_pdb_residue_interactions_distinct_bl_ec["bound_molecule_id"] + "_" + interpro_pdb_residue_interactions_distinct_bl_ec["bound_ligand_id"]
        interpro_pdb_residue_interactions_distinct_bl_ec["type"] = "ligand"

        cath_pdb_residue_interactions_distinct_bl_ec.to_csv(f"{args.outdir}/cath_pdb_residue_interactions_distinct_bl_ec.csv.gz", index = False, compression = "gzip")
        scop_pdb_residue_interactions_distinct_bl_ec.to_csv(f"{args.outdir}/scop_pdb_residue_interactions_distinct_bl_ec.csv.gz", index = False, compression = "gzip")
        interpro_pdb_residue_interactions_distinct_bl_ec.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_bl_ec.csv.gz", index = False, compression = "gzip")
        
        bound_molecules_ligands = pd.concat([cath_pdb_residue_interactions_distinct_bl_ec[["bound_molecule_id", "bound_ligand_id"]].drop_duplicates(), scop_pdb_residue_interactions_distinct_bl_ec[["bound_molecule_id", "bound_ligand_id"]].drop_duplicates(),
            interpro_pdb_residue_interactions_distinct_bl_ec[["bound_molecule_id", "bound_ligand_id"]].drop_duplicates()])
        bound_molecules_ligands.to_csv(f"{args.outdir}/bound_molecules_ligands.csv.gz", index = False, compression = "gzip")
    else:
        print("Loading bound_molecules_ligands")
        bound_molecules_ligands = pd.read_csv(f"{args.outdir}/bound_molecules_ligands.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

    if not os.path.exists(f"{args.outdir}/bound_molceules_sugars.csv.gz"):
        print("Updating sugar EC records")
        cath_pdb_residue_interactions_distinct_sugar_ec = get_updated_enzyme_records(cath_pdb_residue_interactions_distinct_sugar, ec_records_df)
        scop_pdb_residue_interactions_distinct_sugar_ec = get_updated_enzyme_records(scop_pdb_residue_interactions_distinct_sugar, ec_records_df)
        interpro_pdb_residue_interactions_distinct_sugar_ec = get_updated_enzyme_records(interpro_pdb_residue_interactions_distinct_sugar, ec_records_df)

        cath_pdb_residue_interactions_distinct_sugar_ec["uniqueID"] = cath_pdb_residue_interactions_distinct_sugar_ec["bound_molecule_id"] + "_se" + cath_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"]
        cath_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"] = cath_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"].astype(int)
        cath_pdb_residue_interactions_distinct_sugar_ec["type"] = "sugar"
        scop_pdb_residue_interactions_distinct_sugar_ec["uniqueID"] = scop_pdb_residue_interactions_distinct_sugar_ec["bound_molecule_id"] + "_se" + scop_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"]
        scop_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"] = scop_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"].astype(int)
        scop_pdb_residue_interactions_distinct_sugar_ec["type"] = "sugar"
        interpro_pdb_residue_interactions_distinct_sugar_ec["uniqueID"] = interpro_pdb_residue_interactions_distinct_sugar_ec["bound_molecule_id"] + "_se" + interpro_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"]
        interpro_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"] = interpro_pdb_residue_interactions_distinct_sugar_ec["ligand_entity_id_numerical"].astype(int)
        interpro_pdb_residue_interactions_distinct_sugar_ec["type"] = "sugar"

        cath_pdb_residue_interactions_distinct_sugar_ec.to_csv(f"{args.outdir}/cath_pdb_residue_interactions_distinct_sugar_ec.csv.gz", index = False, compression = "gzip")
        scop_pdb_residue_interactions_distinct_sugar_ec.to_csv(f"{args.outdir}/scop_pdb_residue_interactions_distinct_sugar_ec.csv.gz", index = False, compression = "gzip")
        interpro_pdb_residue_interactions_distinct_sugar_ec.to_csv(f"{args.outdir}/interpro_pdb_residue_interactions_distinct_sugar_ec.csv.gz", index = False, compression = "gzip")

        #verifying that the sugar bound ligand ids are not in the bound ligand dataframes.
        assert(len(cath_pdb_residue_interactions_distinct_sugar_ec.loc[cath_pdb_residue_interactions_distinct_sugar_ec.bound_ligand_id.isin(bound_molecules_ligands.bound_ligand_id.unique())]) == 0)
        assert(len(scop_pdb_residue_interactions_distinct_sugar_ec.loc[scop_pdb_residue_interactions_distinct_sugar_ec.bound_ligand_id.isin(bound_molecules_ligands.bound_ligand_id.unique())]) == 0)
        assert(len(interpro_pdb_residue_interactions_distinct_sugar_ec.loc[interpro_pdb_residue_interactions_distinct_sugar_ec.bound_ligand_id.isin(bound_molecules_ligands.bound_ligand_id.unique())]) == 0)

        bound_molecules_sugars = pd.concat([
            cath_pdb_residue_interactions_distinct_sugar_ec[["pdb_id", "bound_molecule_id", "ligand_entity_id", "ligand_entity_description", "ligand_entity_id_numerical", "protein_entity_ec"]].drop_duplicates(), 
            scop_pdb_residue_interactions_distinct_sugar_ec[["pdb_id", "bound_molecule_id", "ligand_entity_id", "ligand_entity_description", "ligand_entity_id_numerical", "protein_entity_ec"]].drop_duplicates(),
            interpro_pdb_residue_interactions_distinct_sugar_ec[["pdb_id", "bound_molecule_id", "ligand_entity_id", "ligand_entity_description", "ligand_entity_id_numerical", "protein_entity_ec"]].drop_duplicates()
            ])
        bound_molecules_sugars.to_csv(f"{args.outdir}/bound_molecules_sugars.csv.gz", index = False, compression = "gzip")
    else:
        print("Loading bound_molecules_sugars")
        bound_molecules_sugars = pd.read_csv(f"{args.outdir}/bound_molecules_sugars.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

    if not os.path.exists(f"{args.outdir}/bound_molecules_sugars_ec_exploded.csv.gz"):
        print("retrieving sugar CIF records")

        cif_ids = bound_molecules_sugars.pdb_id.unique()
        Path(f"{args.outdir}/sugar_cifs").mkdir(parents=True, exist_ok=True)
        for cif_id in cif_ids:
            cif_path = f"{args.outdir}/sugar_cifs/{cif_id}_updated.cif"
            if not os.path.exists(cif_path):
                response = requests.get(f'http://www.ebi.ac.uk/pdbe/entry-files/download/{cif_id}_updated.cif')
                
                with open(cif_path, 'wb') as fp:
                    fp.write(response.content)

        
        branched_entity_list = []
        reader = CifFileReader()
        for cif_id in cif_ids:
            cif_path = f'sugar_cifs/{cif_id}_updated.cif'
            cif_dict = reader.read(cif_path, output='cif_dictionary')
            cif_df = pd.DataFrame(cif_dict)
            if "_pdbx_entity_branch_descriptor" in cif_df.index:
                branch_data = cif_df.loc["_pdbx_entity_branch_descriptor"].values[0]
                if type(next(iter(branch_data.values()))) == str:  # If only one descriptor
                    branched_entity_descriptors = pd.DataFrame({key: [value] for key, value in branch_data.items()})  # Convert scalar value to DataFrame
                else:  # If multiple descriptors
                    branched_entity_descriptors = pd.DataFrame(branch_data)
                branched_entity_descriptors["pdb_id"] = cif_id

                branched_entity_list.append(branched_entity_descriptors)
                
        branched_entity_df = pd.concat(branched_entity_list)
        branched_entity_df["entity_id"] = branched_entity_df.entity_id.astype("int")
        branched_entity_df.loc[branched_entity_df.type == "WURCS"].pdb_id.nunique()

        sugar_wurcs = branched_entity_df.loc[branched_entity_df.type == "WURCS"].groupby(["pdb_id", "entity_id"]).head(1).reset_index()

        sugar_wurcs.rename(columns = {"descriptor" : "WURCS"}, inplace = True)

        bound_molecules_sugars = bound_molecules_sugars.merge(sugar_wurcs[["entity_id", "WURCS", "pdb_id"]], left_on = ["pdb_id", "ligand_entity_id_numerical"], right_on = ["pdb_id", "entity_id"], how = "left", indicator = True)

        bound_molecules_sugars.loc[bound_molecules_sugars._merge != "both", "WURCS"] = "WURCS not available"
        bound_molecules_sugars.drop(columns = ["_merge"], inplace = True)

        bound_molecules_sugars_ec["ec_list"] = bound_molecules_sugars_ec.ec_list.str.split(",")
        bound_molecules_sugars_ec = bound_molecules_sugars_ec.explode("ec_list")
        bound_molecules_sugars_ec.drop(columns = "protein_entity_ec", inplace = True)
        bound_molecules_sugars_ec.to_csv(f"{args.outdir}/bound_molecules_sugars_ec_exploded.csv.gz", compression = "gzip")
    else:
        print("Loading bound_molecules_sugars_ec_exploded")
        bound_molecules_sugars_ec = pd.read_csv(f"{args.outdir}/bound_molecules_sugars_ec_exploded.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

    if not os.path.exists(f"{args.outdir}/bound_molecules_sugars_smiles.pkl"):
        print("retrieving sugar smiles")

        bound_sugars_to_score = bound_molecules_sugars_ec.loc[bound_molecules_sugars_ec.WURCS != "WURCS not available", ["ligand_entity_description","WURCS", "ec_list"]].drop_duplicates()
        bound_sugars_to_score = bound_sugars_to_score.groupby(["ligand_entity_description","WURCS"]).agg({"ec_list": set}).reset_index()

        bound_sugars_to_score["glycoct"] = bound_sugars_to_score["WURCS"].apply(lambda x: get_glycoct_from_wurcs(x))
        bound_sugars_to_score = bound_sugars_to_score.loc[bound_sugars_to_score.glycoct.isna() == False]
        bound_sugars_to_score.to_pickle(f"{args.outdir}/bound_sugars_to_score_temp.pkl") #needed when doing csdb work locally 

        bound_sugars_to_score["csdb"] = bound_sugars_to_score["glycoct"].apply(lambda x: get_csdb_from_glycoct(x))
        bound_sugars_to_score["descriptor"] = bound_sugars_to_score["csdb"].apply(lambda x: get_smiles_from_csdb(x))

        bound_sugars_to_score = pd.read_pickle(f"{args.outdir}/bound_sugars_to_score_temp2.pkl") #need when doing csdb work locally

        bound_sugars_to_score = bound_sugars_to_score.loc[bound_sugars_to_score.descriptor.isna() == False]

        bound_sugars_to_score = bound_sugars_to_score.reset_index()
        bound_sugars_to_score.drop(columns = ["index"], inplace = True)
        bound_sugars_to_score = bound_sugars_to_score.reset_index().rename(columns = {"index": "ligand_index"})

        bound_molecules_sugars_ec = bound_molecules_sugars_ec.merge(bound_sugars_to_score[["ligand_entity_description", "ligand_index", "WURCS", "descriptor"]], on = ["ligand_entity_description","WURCS"], how = "left")

        bound_sugars_to_score["bl_name"] = bound_sugars_to_score["ligand_entity_description"]
        bound_sugars_to_score.rename(columns = {"ligand_index": "ligand_entity_id"}, inplace = True) #do this to run sugars in parity calcs
        bound_sugars_to_score.to_pickle(f"{args.outdir}/bound_sugars_to_score.pkl")

        missing_ligand_index = bound_molecules_sugars_ec.loc[bound_molecules_sugars_ec.descriptor.isna(), ["pdb_id", "entity_id"]].drop_duplicates()
        missing_ligand_index["missing_ligand_index"] = missing_ligand_index.reset_index(drop=True).reset_index().index + bound_molecules_sugars_ec.ligand_index.max() + 1

        bound_molecules_sugars_ec = bound_molecules_sugars_ec.merge(missing_ligand_index, on = ["pdb_id", "entity_id"], how = "left")
        bound_molecules_sugars_ec["ligand_index"].fillna(bound_molecules_sugars_ec["missing_ligand_index"], inplace=True)
        bound_molecules_sugars_ec.drop(columns = "missing_ligand_index", inplace = True)
        bound_molecules_sugars_ec["descriptor"].fillna("SMILES unavailable", inplace = True)

        bound_molecules_sugars_ec.to_pickle(f"{args.outdir}/bound_molecules_sugars_smiles.pkl")
    else:
        print("Loading bound_molecules_sugars_smiles")
        bound_molecules_sugars_ec = pd.read_pickle(f"{args.outdir}/bound_molecules_sugars_smiles.pkl")

    if not os.path.exists(f"{args.outdir}/bound_ligands_to_score.pkl"):
        print("retrieving ligand smiles")
        #should potentially move this part of the script to the yaml fiel that others come from.
        if not os.path.exists(f"{args.outdir}/all_chem_descriptors_bm_ec.csv.gz"):
            if not os.path.exists(f"{args.outdir}/all_chem_descriptors_bm.csv.gz"):
                all_chem_descriptors_query_bm = """
                MATCH (d:ChemicalComponentDesc)<-[g:DESCRIBED_BY]-(cc:ChemicalComponent)<-[h:IS_A]-(p:Entity)<-[i:HAS_ENTITY]-(a:Entry)-[:HAS_ENTITY]->(e:Entity),
                (bl:BoundLigand)<-[j:IS_AN_INSTANCE_OF]-(p)<-[k:HAS_ENTITY]-(bm:BoundMolecule)<-[part:IS_PART_OF]-(bl)
                WHERE p.TYPE = "b" AND p.POLYMER_TYPE = "B" AND
                NOT bl.CHEM_COMP_ID in ["UNX", "UNL"] AND
                e.TYPE = 'p'
                AND e.POLYMER_TYPE = 'P'
                AND e.EC IS NOT NULL
                AND cc.ID = bl.CHEM_COMP_ID
                RETURN DISTINCT 
                    a.ID as pdb_entry_id,
                    bm.UNIQID as bm_id,
                    part.SYM_OPERATOR as bl_sym,
                    bl.UNIQID as bl_id,
                    bl.CHEM_COMP_ID as bl_name,
                    p.UNIQID as ligand_entity_id,
                    p.DESCRIPTION as ligand_entity_description,
                    d.DESCRIPTOR as descriptor, 
                    d.TYPE as descriptor_type, 
                    e.EC as protein_polymer_EC
                """

                all_chem_descriptors_bm = pd.DataFrame([dict(_) for _ in conn.query(all_chem_descriptors_query_bm, db='neo4j')])
                all_chem_descriptors_bm.to_csv(f"{args.outdir}/all_chem_descriptors_bm.csv.gz", compression = "gzip")
            else:
                print("Loading all_chem_descriptors_bm")
                all_chem_descriptors_bm = pd.read_csv(f"{args.outdir}/all_chem_descriptors_bm.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

            all_chem_descriptors_ec = get_updated_enzyme_records(all_chem_descriptors_bm, ec_records_df, ec_col = "protein_polymer_EC")
            all_chem_descriptors_ec["ec_list"] = all_chem_descriptors_ec.ec_list.str.split(",")
            all_chem_descriptors_ec = all_chem_descriptors_ec.explode("ec_list")
            all_chem_descriptors_ec.drop(columns = "protein_polymer_EC", inplace = True)
            all_chem_descriptors_ec.to_csv(f"{args.outdir}/all_chem_descriptors_bm_ec.csv.gz", compression = "gzip")
        else:
            print("Loading all_chem_descriptors_bm_ec")
            all_chem_descriptors_ec = pd.read_csv(f"{args.outdir}/all_chem_descriptors_bm_ec.csv.gz", compression = "gzip", na_values = ["NaN", "None"], keep_default_na = False)

        all_chem_descriptors_smiles = all_chem_descriptors_ec.loc[all_chem_descriptors_ec.descriptor_type == "SMILES_CANONICAL"]
        all_chem_descriptors_smiles_unique_pairs = all_chem_descriptors_smiles.drop_duplicates(["bl_id","ec_list"], keep='first') #get the unique pairs of inchi descriptors and EC numbers

        bound_ligand_descriptors = all_chem_descriptors_smiles_unique_pairs.loc[
            (all_chem_descriptors_smiles_unique_pairs.bl_id.isin(bound_molecules_ligands.bound_ligand_id.unique()))]

        bound_ligands_to_score = bound_ligand_descriptors[["ligand_entity_description", "bl_name", "descriptor", "ec_list"]].drop_duplicates()
        bound_ligands_to_score = bound_ligands_to_score.groupby(["bl_name", "descriptor"]).agg({"ec_list": set, "ligand_entity_description": "first"}).reset_index()
        bound_ligands_to_score = bound_ligands_to_score.reset_index().rename(columns = {"index" : "ligand_entity_id"})
        bound_ligands_to_score["ligand_entity_id"] = bound_ligands_to_score["ligand_entity_id"] + bound_molecules_sugars_ec.ligand_index.max() + 1 #plus one because of 0 index to avoid overlaps
        bound_ligands_to_score.to_pickle(f"{args.outdir}/bound_ligands_to_score.pkl")
    else:
        print("Loading bound_ligands_to_score")
        bound_ligands_to_score = pd.read_pickle(f"{args.outdir}/bound_ligands_to_score.pkl")


    #now do a summary of all the data before exiting.

if __name__ == "__main__":
    main()