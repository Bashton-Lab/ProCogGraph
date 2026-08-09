#!/usr/bin/env python

import re
import numpy as np
import pandas as pd
import requests
from urllib.parse import quote
from bs4 import BeautifulSoup
import json
from Bio.ExPASy import Enzyme as EEnzyme
from pdbeccdutils.helpers.mol_tools import fix_molecule
from rdkit import Chem
from rdkit.Chem import PandasTools
import gzip
import xml.etree.ElementTree as ET
import signal
import glyles
from wurcs_to_iupac import translate as translate_wurcs_to_iupac

#make this function be applicable to extract_pdbe_info script too - need to check if the grouped output at end is appropriate.
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
    ec_records_df_grouped = ec_records_df_grouped.loc[ec_records_df_grouped.DE != "Deleted entry."].copy()
    
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

    return ec_records_df_grouped

#this is a core function - should be stored elsewhere.
#a limitation of this function is that it selects only one representative EC number for transferred entries
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

def get_smiles_from_wurcs_offline(wurcs, timeout_seconds = 15):
    """
    Convert a WURCS glycan descriptor straight to SMILES, entirely offline:
    glypy (WURCS parsing) -> wurcs_to_iupac.translate (IUPAC-condensed
    translation) -> glyles (SMILES generation). No network calls, unlike
    get_glycoct_from_wurcs/get_csdb_from_glycoct/get_smiles_from_csdb above.

    See docs/iupac_translator_plan.md for the full validation history -
    benchmarked at ~89% success on the real cognate-ligand master set and
    ~96% on real PDB-deposited glycans, vs. 0% for the live glycoct-based
    chain above (GlyTouCan's API no longer returns glycoct at all).

    GlyLES is known to hang on malformed input instead of failing fast
    (an ANTLR grammar issue, not specific to this translator), hence the
    timeout - returns np.nan rather than blocking indefinitely.
    """
    if wurcs is np.nan or wurcs is None:
        return np.nan
    try:
        iupac = translate_wurcs_to_iupac(wurcs)
    except Exception:
        return np.nan

    def _timeout_handler(signum, frame):
        raise TimeoutError("glyles.convert timed out")

    previous_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(timeout_seconds)
    try:
        result = glyles.convert(glycan = iupac, verbose = None)
    except Exception:
        return np.nan
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)

    if not result or not result[0][1]:
        return np.nan
    return result[0][1]

def pdbe_sanitise_smiles(smiles, return_mol = False, return_sanitisation = False):
    """
    Sanitises a smiles string using pdbeccdutils fix_molecule functions and
    returns a canonical smiles string or RDKit molecule object. Requires that
    smiles string can be loaded into an RDKit molecule object - returns none if
    this is not possible.
    """
    try:
        mol = Chem.MolFromSmiles(Chem.CanonSmiles(smiles))
    except:
        if return_sanitisation:
            return None, None
        else:
            return None
    
    sanitisation = fix_molecule(mol)
    if sanitisation:
        if return_mol:
            if return_sanitisation:
                return mol, sanitisation
            else:
                return mol
        else:
            sanitised_smiles = Chem.CanonSmiles(Chem.MolToSmiles(mol))
            if return_sanitisation:
                return sanitised_smiles, sanitisation
            else:
                return sanitised_smiles
    else:
        sanitised_smiles = None
        if return_sanitisation:
            return sanitised_smiles, sanitisation
        else:
            return sanitised_smiles
        
def extract_interpro_domain_annotations(xml_file):
    with gzip.open(xml_file, 'rb') as f:
        tree = ET.parse(f)
        root = tree.getroot()

        # Find all interpro elements in the XML file
        interpro_elements = root.findall('.//interpro')

        interpro_info = []
        # Iterate through each interpro element
        for interpro in interpro_elements:
            superfamily_annotations = {}
            gene3d_annotations = {}
            interpro_id = interpro.attrib['id']
            interpro_short_name = interpro.attrib['short_name']
            # Store PFAM annotations for the interpro ID
            interpro_info.append({"interpro_accession": interpro_id,
                                        "interpro_name": interpro_short_name})
    interpro_info_df = pd.DataFrame(interpro_info)       

    return interpro_info_df

def get_scop_domains_info(domain_info_file, descriptions_file):
    
    def clean_and_merge_scop_col(df, column_id, description_df):
        level = df[column_id].str.split("=").str.get(0).values[0]
        df[column_id] = df[column_id].str.split("=").str.get(1).astype(int)
        df = df.merge(description_df.loc[description_df.level == level, ["level_sunid", "level", "level_description"]],left_on = column_id, right_on = "level_sunid", indicator = True)
        df.rename(columns = {"level_description": f"{level}_description"}, inplace = True)
        assert len(df.loc[df._merge != "both"]) == 0
        df.drop(columns = ["_merge", "level_sunid", "level"], inplace = True)
        return df
    
    scop_domains_info = pd.read_csv(domain_info_file, sep = "\t", comment = "#", header = None, names = ["scop_id", "pdb_id", "scop_description", "sccs", "domain_sunid", "ancestor_sunid"])
    scop_id_levels = ["cl_id", "cf_id", "sf_id", "fa_id", "dm_id", "sp_id", "px_id"]
    scop_domains_info[scop_id_levels] = scop_domains_info.ancestor_sunid.str.split(",", expand = True)
    scop_descriptions = pd.read_csv(descriptions_file, sep = "\t", comment = "#" , header = None, names = ["level_sunid", "level", "level_sccs", "level_sid", "level_description"])

    for column in scop_id_levels:
        scop_domains_info = clean_and_merge_scop_col(scop_domains_info, column, scop_descriptions)
    
    scop_domains_info.drop(columns = ["pdb_id", "scop_description"], inplace = True)
    return scop_domains_info

def get_pfam_annotations(pfam_clans_file):
    """Pfam consolidated the three files this used to read separately -
    pfamA.txt.gz (pfam accession/name/description), clan_membership.txt.gz
    (clan-to-pfam membership) and clan.txt.gz (clan details) - into a
    single Pfam-A.clans.tsv.gz (2026 restructure, no header row):
    pfam_accession, clan_acc, clan_id, pfam_name, pfam_description. Pfam
    families with no clan have empty clan_acc/clan_id fields (parsed as
    NaN by pandas). See docs/reference_data_download_plan.md #2.

    Schema note: the consolidated file has no equivalent of the old
    clan.txt.gz's free-text `clan_comment` field - Pfam no longer
    publishes a bulk clan comment/description beyond the short clan id
    (e.g. "GPCR_A"). clan_description below is that short id;
    clan_comment is left null rather than duplicating clan_description
    into it, since that would misrepresent it as distinct information.
    """
    pfam_clans = pd.read_csv(pfam_clans_file, sep = "\t", header = None,
        names = ["pfam_accession", "clan_acc", "clan_id", "pfam_name", "pfam_description"])
    pfam_clans["clan"] = pfam_clans["clan_acc"]
    pfam_clans["pfam"] = pfam_clans["pfam_accession"].where(pfam_clans["clan_acc"].notna())
    pfam_clans["clan_description"] = pfam_clans["clan_id"]
    pfam_clans["clan_comment"] = np.nan
    return pfam_clans[["pfam_accession", "pfam_name", "pfam_description", "clan", "pfam", "clan_acc", "clan_description", "clan_comment"]]

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

def get_chem_comp_descriptors(ccd_doc, comp_id_list):
    """Resolve a list of PDB chemical component codes to a single SMILES
    descriptor each, using the CCD's own _pdbx_chem_comp_descriptor loop
    (OpenEye descriptors preferred when present, else the first
    RDKit-parseable SMILES row). Moved here from process_all_pdb_contacts.py
    so it can be reused by preprocess_cofactors.py without a cross-script
    import - it's a generic CCD-parsing utility, not specific to the
    contacts pipeline."""
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

def classify_ec_completeness(ec):
    """Returns how many of an EC number's 4 dot-separated segments are
    fully resolved (non "-") before the first wildcard, reading
    left-to-right and stopping at the first "-". E.g. "1.1.1.10" -> 4,
    "1.1.1.-" -> 3, "1.1.-.-" -> 2, "1.-.-.-" -> 1. Used by
    preprocess_cofactors.py to decide whether a source's EC value is exact,
    safely broadcastable (subsubclass-level, 3), or too coarse to use at
    all (class/subclass-level, 1 or 2) - see docs/cofactor_coverage_plan.md
    ("The broadcast problem, and its fix") for why levels 1 and 2 are
    excluded: broadcasting them down to every terminal EC in that class/
    subclass would claim structurally unrelated enzymes share a cofactor
    just because they share a leading EC digit or two."""
    segments = ec.split(".")
    level = 0
    for segment in segments:
        if segment == "-":
            break
        level += 1
    return level

def broadcast_subsubclass_ec(ec, terminal_ec_list):
    """Expands a subsubclass-level partial EC (e.g. "1.1.1.-") to every
    matching terminal EC in terminal_ec_list. Reuses the existing
    return_partial_EC_list matching logic (already used elsewhere in this
    pipeline for SIFTS' partial EC annotations) rather than re-implementing
    prefix matching. Callers MUST have already checked
    classify_ec_completeness(ec) == 3 before calling this - it will happily
    (and, per the cofactor plan, incorrectly) broadcast a class- or
    subclass-level wildcard too, since return_partial_EC_list itself has no
    concept of "too coarse to use"."""
    matches = return_partial_EC_list(ec, terminal_ec_list)
    return matches if isinstance(matches, list) else []

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
    residue_ec_records_grouped["TRANSFER"] = residue_ec_records_grouped["TRANSFER"].str.strip(",")
    residue_ec_records_grouped.rename(columns = {"TRANSFER" : "ec_list"}, inplace = True)
    
    df_merged = df.merge(residue_ec_records_grouped, on = ec_col, how = "left", indicator = True)
    assert(len(df_merged.loc[df_merged["_merge"] != "both"]) == 0)
    df_merged.drop(columns = "_merge", inplace = True)
    df_merged = df_merged.loc[df_merged.ec_list != ""] #remove any rows where the ec_list is empty - we cant process these anyway.
    return(df_merged)

def parse_cddf(file_path, domain_list):
    with open(file_path, 'r') as file:
        entries = {}
        lines = []
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                continue
            if line.startswith("//"):
                entries[domain] = lines
                lines = []
            elif line.startswith("DOMAIN"):
                domain = line[10:]
                lines.append(line)
            else:
                lines.append(line)
    matching_entries = []
    for domain in domain_list:
        if domain in entries.keys():
            current_entry = {}
            current_segment = {}
            for line in entries[domain]:
                if line.startswith('SEGMENT'):
                    if 'SEGMENTS' not in current_entry:
                        current_entry['SEGMENTS'] = []
                    if current_segment:
                        current_entry['SEGMENTS'].append(current_segment)
                    current_segment = {'SEGMENT': line.split()[1]}
                elif line.startswith('ENDSEG'):
                    if current_segment:
                        current_entry['SEGMENTS'].append(current_segment)
                    current_segment = {}
                else:
                    tag = line[0:9].strip()
                    value = line[10:]

                    if tag == 'FORMAT':
                        if 'FORMAT' not in current_entry:
                            current_entry['FORMAT'] = []
                        current_entry['FORMAT'].append(value)
                    elif tag in ['DOMAIN', 'PDBID', 'VERSION', 'VERDATE', 'NAME', 'SOURCE', 'CATHCODE', 'CLASS', 'ARCH', 'TOPOL', 'HOMOL', 'DLENGTH', 'DSEQH', 'DSEQS', 'NSEGMENTS']:
                        if tag not in current_entry:
                            current_entry[tag] = []
                        current_entry[tag].append(value)
                    elif tag in ['SRANGE', 'SLENGTH', 'SSEQH', 'SSEQS']:
                        current_segment[tag] = value
                    elif tag == 'COMMENTS':
                        if 'COMMENTS' not in current_entry:
                            current_entry['COMMENTS'] = []
                        current_entry['COMMENTS'].append(value)
                    else:
                        # Invalid tag
                        continue
            matching_entries.append(current_entry)
    return matching_entries


def build_cath_dataframe(parsed_data):
    combine_keys = ["DOMAIN", "PDBID", "VERSION", "VERDATE", "NAME", "SOURCE", "CATHCODE", "CLASS", "ARCH", "TOPOL", "HOMOL", "DLENGTH", "DSEQH", "DSEQS", "NSEGMENTS"]
    cath_df_columns = {"DOMAIN": "cath_domain" , "VERSION": "cath_db_version", "VERDATE": "cath_db_verdate", "NAME":"cath_name", "SOURCE": "cath_source", "CATHCODE": "cath_code", "CLASS": "cath_class_name", "ARCH":"cath_architecture_name", "TOPOL": "cath_topology_name", "HOMOL": "cath_homologous_superfamily_name", "DLENGTH": "cath_domain_length", "DSEQH":"cath_domain_seq_header", "DSEQS": "cath_domain_seqs", "NSEGMENTS": "cath_num_segments", "SEGMENTS": "cath_segments_dict"}
    dfs = []
    for entry in parsed_data:
        df_dict = {}
        for key, values in entry.items():
            if key in combine_keys:
                df_dict[cath_df_columns[key]] = ' '.join(values) if isinstance(values, list) else values
            elif key == "FORMAT":
                continue
            else:
                df_dict[cath_df_columns[key]] = values
        dfs.append(pd.DataFrame([df_dict]))
    dfs_combined = pd.concat(dfs, ignore_index=True)
    topology_regex = r'^(\d+\.\d+\.\d+)\.'
    architecture_regex = r'^(\d+\.\d+)\.'
    class_regex = r'^(\d+)\.'
    
    dfs_combined["cath_class"] = dfs_combined.cath_code.str.extract(class_regex, expand = True)
    dfs_combined["cath_architecture"] = dfs_combined.cath_code.str.extract(architecture_regex, expand = True)
    dfs_combined["cath_topology"] = dfs_combined.cath_code.str.extract(topology_regex, expand = True)
    dfs_combined["cath_homologous_superfamily"] = dfs_combined["cath_code"]
    return dfs_combined

def build_g3dsa_dataframe(cath_names, cath_domain_list): #unlike cath, g3dsa annotations are at superfamily level (also works for xml cath annotations from sifts which are at homolsuperfam level)
    topology_regex = r'^(\d+\.\d+\.\d+)\.'
    architecture_regex = r'^(\d+\.\d+)\.'
    class_regex = r'^(\d+)\.'
    domain_list = []
    for homologous_superfamily in cath_domain_list:
        homologous_superfamily_name = cath_names.loc[cath_names.cath_code == homologous_superfamily, "name"].values[0]
        topology = re.search(topology_regex, homologous_superfamily).group(1)
        topology_name = cath_names.loc[cath_names.cath_code == topology, "name"].values[0]
        architecture = re.search(architecture_regex, homologous_superfamily).group(1)
        architecture_name = cath_names.loc[cath_names.cath_code == architecture, "name"].values[0]
        class_ = re.search(class_regex, homologous_superfamily).group(1)
        class_name = cath_names.loc[cath_names.cath_code == class_, "name"].values[0]
        domain = {"cath_domain": homologous_superfamily, "cath_name" : homologous_superfamily_name, "cath_code": homologous_superfamily, "cath_homologous_superfamily": homologous_superfamily, "cath_homologous_superfamily_name": homologous_superfamily_name, "cath_topology": topology, "cath_topology_name": topology_name,  "cath_architecture": architecture, "cath_architecture_name": architecture_name, "cath_class": class_, "cath_class_name" : class_name}
        domain_list.append(domain)
    domain_df = pd.DataFrame(domain_list)
    return domain_df

def get_scop2_domains_info(domain_info_file, descriptions_file):
    scop_domains_info = pd.read_csv(domain_info_file, sep = " ", comment = "#", header = None, names = ["FA-DOMID", "FA-PDBID","FA-PDBREG","FA-UNIID","FA-UNIREG","SF-DOMID","SF-PDBID","SF-PDBREG","SF-UNIID","SF-UNIREG","SCOPCLA"])
    scop_fa_domains_info = scop_domains_info[["FA-DOMID", "SCOPCLA"]].groupby("FA-DOMID").agg({"SCOPCLA": list}).reset_index()
    scop_fa_domains_info["SCOPCLA"] = scop_fa_domains_info["SCOPCLA"].str.join(";")

    scop_sf_domains_info = scop_domains_info[["SF-DOMID", "SCOPCLA"]]
    scop_sf_domains_info["SCOPCLA"] = scop_sf_domains_info["SCOPCLA"].str.extract("(.*),FA=.*$", expand = True) #remove the family level from sueprfamily level domains
    scop_sf_domains_info = scop_sf_domains_info.groupby("SF-DOMID").agg({"SCOPCLA": list}).reset_index()
    scop_sf_domains_info["SCOPCLA"] = scop_sf_domains_info["SCOPCLA"].str.join(";")
    with open(descriptions_file) as file:
        rows = []
        for row in file:
            if not row.startswith("#"):
                pattern = re.compile(r"(\d+) (.+)")
                matches = pattern.match(row)
                rows.append([matches.group(1), matches.group(2)])
    scop_descriptions = pd.DataFrame(rows, columns = ["NODE_ID", "NODE_NAME"])

    return scop_sf_domains_info, scop_fa_domains_info, scop_descriptions