#!/usr/bin/env python

import argparse
from pathlib import Path
import scipy.stats as stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



def parse_args():
    parser = argparse.ArgumentParser("Compute lambda kistler")
    parser.add_argument("length", type=Path, help="Path to mapdamage lengthdistrib file")
    parser.add_argument("output", type=str, help="Output file basename")
    parser.add_argument("--read_len_min", type=int, default=30, dest = 'read_len', help="Minimum read length to consider")

    return parser.parse_args()

def compute_kistler_lambda(length, min_len, basename):
    df = pd.read_table(length, comment='#')
    df.columns = ['Std', 'Length', 'Occurences']
    read_lengths = np.repeat(df['Length'], df['Occurences'])

    # Find the mode (peak of the distribution)
    hist_values, bin_edges = np.histogram(read_lengths, bins=50, density=True)
    mode_index = np.argmax(hist_values)  # Index of peak
    mode_value = bin_edges[mode_index]  # Mode (starting point for decreasing part)

    # Select only the decreasing part of the distribution
    filtered_data = read_lengths[read_lengths > max(min_len, mode_value)]

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

    return mode_value, lambda_hat

def geom_pmf(k, p):
    return stats.geom.pmf(k, p)

def compute_geom_param(length, read_len_min, basename):
    df = pd.read_table(length, comment="#")
    df.columns = ['Std', 'Length', 'Occurences']
    df = df.groupby('Length')['Occurences'].sum().reset_index().sort_values('Length')

    fits = dict()
    for i in range(read_len_min, read_len_min + 5):
        _ = df.query('Length >= @i')
        x = np.array(_.Length) - i
        y = np.array(_.Occurences)/_.Occurences.sum()
        d = dict()

        d['popt'], d['pcov'], d['infodict'], d['mesg'], d['ier'] = curve_fit(geom_pmf, x, y, bounds=((0),(1)), full_output=True)
        fits[i] = d
    best_fit = min(fits, key=lambda x: np.absolute(fits[x]['infodict']['fvec'].sum()))

    fig = plt.figure()
    plt.plot(list(fits.keys()), [np.absolute(fits[x]['infodict']['fvec'].sum()) for x in fits.keys()])
    plt.suptitle(f"{basename} Sum of residuals at each minimum read length")
    plt.title(f"Best fit: {best_fit} bp")
    fig.savefig(f"{basename}_geom_mode_best_fit.png")

    fig2 = plt.figure()
    plt.plot(x+best_fit, y, 'o', label='observed')
    fit_y = geom_pmf(x, fits[best_fit]['popt'][0])
    fit_err = np.sqrt(np.diag(fits[best_fit]['pcov']))
    plt.plot(x+best_fit,fit_y , 'r-', label = f"fitted geom pmf(p={round(fits[best_fit]['popt'][0], 3)}, loc={best_fit})")
    plt.fill_between(x+best_fit, fit_y - fit_err, fit_y + fit_err, color='r', alpha=0.5, label='+/- 1 std')
    plt.legend()
    plt.suptitle(f"{basename} Geometric PMF fit")
    fig2.savefig(f"{basename}_geom_pmf_fit.png")

    return best_fit, fits[best_fit]['popt'][0]





if __name__ == "__main__":
    args = parse_args()
    mode_expo, la = compute_kistler_lambda(args.length, args.read_len, args.output)
    mode_geom, p = compute_geom_param(args.length, args.read_len, args.output)

    with open(f"{args.output}_params.tsv", 'w') as f:
        f.write("mode_expo\tlambda\tmode_geom\tp\n")
        f.write(f"{mode_expo}\t{la}\t{mode_geom}\t{p}\n")
