import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'TO DO')
parser.add_argument('--chunk_list', metavar = '', nargs='+', 
    help = "")
parser.add_argument('--ligands_score_subset', metavar = '', type = str,
                    help = "")
parser.add_argument('--output_pkl', metavar = '', type = str,
    help = "")

args = parser.parse_args()

all_chunks = args.chunk_list
all_chunks_df = pd.DataFrame()
for chunk in all_chunks:
    chunk_df = pd.DataFrame(pd.read_pickle(chunk))
    all_chunks_df = pd.concat([all_chunks_df,chunk_df])

all_chunks_df = all_chunks_df.loc[all_chunks_df.error.isna()]
all_chunks_df["ec"] = all_chunks_df["ec"].str.split(",")
all_chunks_df.to_pickle(args.output_pkl)