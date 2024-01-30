#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser("Create preseq plot file")
    parser.add_argument("preseq_curve", type=Path, help="Path to preseq curve file")
    parser.add_argument("-o", type=Path, dest="output", help="Output file")
    parser.add_argument("-s", type=str, dest="sample", help="Sample name")
    parser.add_argument("-p", type=str, dest="preseq_mode", help="Preseq mode")
    return parser.parse_args()


def read_preseq_curve(preseq_curve):
    df = pd.read_csv(preseq_curve, sep="\t")
    return df


def plot_preseq(psc, sample, preseq_mode, out):
    if preseq_mode == "c_curve":
        x_col = "total_reads"
        y_col = "distinct_reads"
    elif preseq_mode == "lc_extrap":
        x_col = "TOTAL_READS"
        y_col = "EXPECTED_DISTINCT"

    plt.figure(figsize=(8, 6))
    plt.plot(
        psc[x_col],
        psc[y_col],
        label="Reads",
        linestyle="--",
        marker="o",
    )
    plt.plot(
        df[x_col],
        df[x_col],
        label="Maximal complexity",
        linestyle="-",
    )
    if preseq_mode == "lc_extrap":
        plt.fill_between(
            df["TOTAL_READS"],
            df["LOWER_0.95CI"],
            df["UPPER_0.95CI"],
            color="skyblue",
            alpha=0.4,
            label="95% Confidence Interval",
        )
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f"Preseq {preseq_mode} curve - {sample}")
    plt.legend()
    plt.grid(True)
    plt.savefig(out)


if __name__ == "__main__":
    args = parse_args()
    df = read_preseq_curve(args.preseq_curve)
    plot_preseq(df, args.sample, args.preseq_mode, args.output)
