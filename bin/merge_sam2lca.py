#!/usr/bin/env python


import argparse
from pathlib import Path
import pandas as pd
from functools import reduce



def parse_args():
    parser = argparse.ArgumentParser(
        "Create sam2lca json file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
            "-r",
            type=str,
            dest="rank",
            default='species',
            help="Taxonomic rank used for merging."
        )

    return parser.parse_args()

def filter_df(df, rank, sample_name):
    df = pd.read_csv(df)
    df = df[df["rank"] == rank][['name','count_descendant']]
    df = df.rename(columns={'name': rank, 'count_descendant': sample_name})
    return df

if __name__ == "__main__":

    args = parse_args()
    cwd = Path.cwd()
    sam2lca_results = cwd.glob("*.sorted.sam2lca.csv")
    dfs = []
    for res in sam2lca_results:
        samp_name = res.stem.split(".")[0]
        dfs.append(filter_df(res, args.rank, samp_name))

    df_merged = reduce(lambda left,right: pd.merge(left,right,on=args.rank, how='outer'), dfs)
    df_merged.fillna(0).to_csv("merged.sam2lca.csv", index=False)


