import argparse
import pandas as pd
import pickle 

parser = argparse.ArgumentParser(description = 'TO DO')
parser.add_argument('--chunk_list', metavar = '', nargs='+', 
    help = "")
parser.add_argument('--ligands_score_subset', metavar = '', type = str,
                    help = "")
parser.add_argument('--all_ligands_pkl', metavar = '', type = str,
                    help = "")
parser.add_argument('--output_pkl', metavar = '', type = str,
    help = "")

args = parser.parse_args()

all_chunks = args.chunk_list
all_chunks_loaded = []
for chunk in all_chunks:
    chunk_df = pd.DataFrame(pd.read_pickle(chunk))
    all_chunks_loaded.append(chunk_df)

all_chunks_df = pd.concat(all_chunks_loaded)

ligands_score_subset = pd.read_pickle(args.ligands_score_subset)
all_ligands = pd.read_pickle(args.all_ligands_pkl)
all_chunks_df = pd.merge(ligands_score_subset, all_chunks_df, left_on = "ligand_entity_id", right_on = "pdb_ligand", indicator = True)
assert(len(all_chunks_df.loc[all_chunks_df._merge != "both"]) == 0)
all_chunks_df.drop(columns = ["ligand_entity_id", "pdb_ligand", "_merge"], inplace = True)
all_chunks_df = pd.merge(all_chunks_df, all_ligands, how = "right", on = ["descriptor", "protein_polymer_EC"], indicator = True)
assert(len(all_chunks_df.loc[all_chunks_df._merge != "both"]) == 0)
all_chunks_df.drop(columns = ["_merge"], inplace = True)

max_indexes = all_chunks_df.groupby('bl_id')['score'].idxmax()

# Mark rows as Cognate based on conditions
all_chunks_df['isCognate'] = False  # Initialize all rows as non-cognate
cognate_mask = (all_chunks_df['score'] > 0.55) & all_chunks_df.index.isin(max_indexes)
all_chunks_df.loc[cognate_mask, 'isCognate'] = True

with open(args.output_pkl, "wb") as file:
    pickle.dump(all_chunks_df, file)