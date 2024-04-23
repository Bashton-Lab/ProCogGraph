


ec_records_df_grouped = process_ec_records(args.enzyme_dat_file , args.enzyme_class_file)
ec_records_df = ec_records_df_grouped.explode("ID")
sifts_chains_ec = sifts_chains.loc[sifts_chains.EC_NUMBER != "?"]
sifts_chains_uniprot = sifts_chains.groupby(["PDB", "CHAIN"]).agg({"ACCESSION": set}).reset_index()
sifts_chains_uniprot["ACCESSION"] = sifts_chains_uniprot["ACCESSION"].apply(lambda x: "|".join(x)) #join the list of uniprot accessions with a pipe for downstream neo4j integration
sifts_chains_uniprot.rename(columns = {"ACCESSION" : "uniprot_accession"}, inplace = True)

sifts_chains_ec = sifts_chains_ec.groupby(["PDB"]).agg({"EC_NUMBER": set}).reset_index() #group these into a set of pdb associated ec's
sifts_chains_ec["EC_NUMBER"] = sifts_chains_ec["EC_NUMBER"].apply(lambda x: ",".join(x)) #join the list of EC numbers into a single string for ec records function 
sifts_chains_ec = get_updated_enzyme_records(sifts_chains_ec, ec_records_df, ec_col = "EC_NUMBER")
sifts_chains_ec.rename(columns = {"EC_NUMBER": "protein_entity_ec", "ACCESSION" : "uniprot_accession"}, inplace = True)

all_pdbs = sifts_chains_ec.PDB.unique()