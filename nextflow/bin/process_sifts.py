#/usr/bin/env python

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sifts_file', type=str, help='')
    parser.add_argument('--assemblies_file', type=str, help='')
    args = parser.parse_args()

    sifts_file = pd.read_csv(args.sifts_file, sep = "\t", skiprows = 1)
    assemblies_file = pd.read_csv(args.assemblies_file, usecols = ["ASSEMBLIES","PREFERED_ASSEMBLIES"]) #from https://ftp.ebi.ac.uk/pub/databases/pdbe-kb/complexes/, see https://www.biorxiv.org/content/10.1101/2023.05.15.540692v1.full
    assemblies_file[["PDB", "ASSEMBLY_ID"]] = assemblies_file["ASSEMBLIES"].str.split("_", expand = True)
    assemblies_file_prefered = assemblies_file.loc[assemblies_file.PREFERED_ASSEMBLIES == 1] #bool true
    sifts_file_ec = sifts_file.loc[~sifts_file.EC_NUMBER.isin(["?"]), ["PDB"]].drop_duplicates()
    sifts_assemblies = sifts_file_ec.merge(assemblies_file_prefered, how = "left", on = "PDB", indicator = True)
    assert(len(sifts_assemblies.loc[sifts_assemblies._merge != "both"]) == 0)
    sifts_assemblies.to_csv("sifts_assemblies.csv", index = False)

if __name__ == "__main__":
    main()