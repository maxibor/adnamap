#!/usr/bin/env python

import argparse
from pathlib import Path
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser("Compute lambda kistler")
    parser.add_argument("length", type=Path, help="Path to mapdamage lengthdistrib file")
    parser.add_argument("output", type=str, help="Output file basename")

    return parser.parse_args()


def compute_lambda(length, basename):
    df = pd.read_table(length, comment='#')
    df.columns = ['Std', 'Length', 'Occurences']
    read_lengths = np.repeat(df['Length'], df['Occurences'])

    # Find the mode (peak of the distribution)
    hist_values, bin_edges = np.histogram(read_lengths, bins=50, density=True)
    mode_index = np.argmax(hist_values)  # Index of peak
    mode_value = bin_edges[mode_index]  # Mode (starting point for decreasing part)

    # Select only the decreasing part of the distribution
    filtered_data = read_lengths[read_lengths > mode_value]

    # Fit an exponential distribution to the decreasing part
    loc, scale = stats.expon.fit(filtered_data, floc=mode_value)  # Force location to 0
    lambda_hat = 1 / scale  # Estimated λ (rate parameter)

    # Generate x values for plotting
    x = np.linspace(filtered_data.min(), filtered_data.max(), 1000)
    pdf_fitted = stats.expon.pdf(x, loc, scale)

    # Plot histogram and fitted exponential distribution
    plt.figure(figsize=(8, 6))
    plt.hist(read_lengths, bins=50, density=True, alpha=0.6, color='black')
    plt.axvline(mode_value, color='blue', linestyle="--", label=f"Mode: {mode_value:.2f}")
    plt.plot(x, pdf_fitted, 'r-', lw=2, label=f'Fitted Exponential (λ={lambda_hat:.4f})')

    plt.xlabel("Read Length")
    plt.ylabel("Density")
    plt.title(f"{basename} read length distribution")
    plt.legend()
    plt.savefig(f"{basename}.png")

    with open(f"{basename}_lambda.tsv", 'w') as f:
        f.write("mode\tlambda\n")
        f.write(f"{mode_value}\t{lambda_hat}\n")

    return mode_value, lambda_hat

if __name__ == "__main__":
    args = parse_args()
    compute_lambda(args.length, args.output)
