from rdfreader import RDFParser
from rdfreader.chem.mol import Molecule
import os
import re   
import argparse
import pandas as pd
from rdkit.Chem import PandasTools
from pathlib import Path
#example usage: python preprocess_rhea.py --rhea_ec_mapping rhea_ec_mapping.tsv --rhea_reaction_directions rhea_reaction_directions.tsv --rd_dir /path/to/rhea/rd/files --outdir /path/to/output/directory --chebi_names /path/to/chebi_names.tsv.gz

def main():
    parser = argparse.ArgumentParser(description='Preprocess Rhea data')
    parser.add_argument('--rhea_ec_mapping', metavar = '', type = str,
        help = "File containing Rhea to EC number mapping")
    parser.add_argument('--rhea_reaction_directions', metavar = '', type = str,
        help = "File containing Rhea reaction directions")
    parser.add_argument('--rd_dir', metavar = '', type = str,
        help = "Directory containing Rhea reaction data files (download from RHEA website)")
    parser.add_argument('--chebi_names', metavar = '', type = str,
        help = "File containing ChEBI name mappings (download from ChEBI website) - expected to be in tsv.gz format")
    parser.add_argument('--outdir', metavar = '', type = str,
        help = "Output directory")

    args = parser.parse_args()
    Path(args.outdir).mkdir(parents = True, exist_ok = True)
    if not os.path.exists(f"{args.outdir}/rhea_rd_parsed.pkl"):
        mol_dict = {}
        unique_id = 0

        for file in os.listdir(f"{args.rd_dir}"):
            if file.endswith(".rd"):
                with open(f"{args.rd_dir}/{file}", "r") as rdf_file:

                    # create a RDFParser object, this is a generator that yields Reaction objects
                    rdfreader = RDFParser(
                        rdf_file,
                        except_on_invalid_molecule=False,  # will return None instead of raising an exception if a molecule is invalid
                        except_on_invalid_reaction=False,  # will return None instead of raising an exception if a reaction is invalid 
                    )
                    for rxn in rdfreader:
                        if rxn is None:
                            continue # the parser failed to read the reaction, go to the next one
                        reaction_id = rxn.id
                        # rxn is a Reaction object, it is several attributes, including:
                        reaction_smiles = rxn.smiles # reaction SMILES string
                        reaction_properties = rxn.properties # a dictionary of properties extracted from the RXN record

                        reactants = rxn.reactants # a list of Molecule objects
                        products = rxn.products


                        for reactant in reactants:
                            reactant_smiles = reactant.smiles
                            reactant_id = reactant.metadata["molecule_name"]
                            mol_dict[unique_id] = {"reaction_id": reaction_id,
                                            "reaction_properties": reaction_properties,
                                            "reaction_smiles": reaction_smiles,
                                            "mol_type": "reactant",
                                            "compound_id": reactant_id,
                                            "smiles": reactant_smiles}
                            unique_id += 1

                        for product in products:
                            product_smiles = product.smiles
                            product_id = product.metadata["molecule_name"]
                            mol_dict[unique_id] = {"reaction_id": reaction_id,
                                            "reaction_properties": reaction_properties,
                                            "reaction_smiles": reaction_smiles,
                                            "mol_type": "product",
                                            "compound_id" : product_id,
                                            "smiles": product_smiles}
                            unique_id += 1

        reactions_df = pd.DataFrame(mol_dict).T
        reactions_df.to_pickle(f"{args.outdir}/rhea_rd_parsed.pkl")
    else:
        reactions_df = pd.read_pickle(f"{args.outdir}/rhea_rd_parsed.pkl")
        
    rhea2ec = pd.read_csv(f"{args.rhea_ec_mapping}", sep = "\t")
    rhea_dir = pd.read_csv(f"{args.rhea_reaction_directions}", sep = "\t")

    reactions_df["reaction_id"] = reactions_df.reaction_id.astype("int")
    rhea2ec["RHEA_ID"] = rhea2ec["RHEA_ID"].astype("int")

    rheamerge = rhea2ec.merge(rhea_dir, left_on = "MASTER_ID" , right_on = "RHEA_ID_MASTER", how = "left")


    reactions_df_merged = reactions_df.merge(rheamerge[["RHEA_ID_LR", "ID"]], left_on = "reaction_id", right_on = "RHEA_ID_LR", how = "inner")
    reactions_df_merged.loc[reactions_df_merged.compound_id.str.startswith("CHEBI"), "COMPOUND_ID"] = reactions_df_merged.loc[reactions_df_merged.compound_id.str.startswith("CHEBI"), "compound_id"].apply(lambda x: re.findall(r"CHEBI:(\d+)", x)[0]) #in chebi names format
    reactions_df_merged.loc[reactions_df_merged.compound_id.str.startswith("CHEBI") == False, "COMPOUND_ID"] = -1
    reactions_df_merged["COMPOUND_ID"] = reactions_df_merged["COMPOUND_ID"].astype("int")

    chebi_names = pd.read_csv(f"{args.chebi_names}", sep = "\t", compression = "gzip")
    #get the first name for each compound ID in the chebi names file
    chebi_names = chebi_names.groupby("COMPOUND_ID").agg({"NAME": "first"}).reset_index()
    print(chebi_names)
    reactions_df_merged = reactions_df_merged.merge(chebi_names, left_on = "COMPOUND_ID", right_on = "COMPOUND_ID", how = "left")
    print(reactions_df_merged)
    reactions_df_merged["NAME"] = reactions_df_merged["NAME"].fillna(reactions_df_merged["compound_id"])
    print(reactions_df_merged)
    reactions_df_merged_filtered = reactions_df_merged.loc[reactions_df_merged.smiles.isna() == False].copy()
    reactions_df_merged_filtered["reaction_id"] = "RHEA:" + reactions_df_merged_filtered["reaction_id"].astype("str")
    reactions_df_merged_filtered.rename(columns = {"ID": "entry", "NAME": "compound_name"}, inplace = True)
    reactions_df_merged_filtered["ligand_db"] = reactions_df_merged_filtered["compound_id"].astype("str")
    reactions_df_merged_filtered.to_pickle(f"{args.outdir}/rhea_reactions_all.pkl")
    reactions_df_merged_filtered_grouped = reactions_df_merged_filtered[['reaction_id', 'smiles', 'entry', 'compound_name', 'ligand_db', 'compound_id']].groupby(["entry", "compound_name", "smiles", "ligand_db", "compound_id"]).agg({"reaction_id": set}).reset_index()
    reactions_df_merged_filtered_grouped["reaction_id"] = reactions_df_merged_filtered_grouped["reaction_id"].str.join("|")
    PandasTools.AddMoleculeColumnToFrame(reactions_df_merged_filtered_grouped, smilesCol='smiles')
    reactions_df_merged_filtered_grouped.to_pickle(f"{args.outdir}/rhea_reactions.pkl")


if __name__ == "__main__":
    main()