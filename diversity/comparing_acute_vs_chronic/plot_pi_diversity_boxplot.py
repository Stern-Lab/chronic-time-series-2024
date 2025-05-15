#!/usr/bin/env python3
"""
Plot comparison of SARS-CoV-2 nucleotide diversity (π) in acute vs. chronic samples.

This script:
- Loads diversity summary CSVs for acute and chronic samples
- Computes a Mann–Whitney U test on replicate-concordant nucleotide diversity
- Produces a log-scaled boxplot with statistical annotation

Usage:
    python plot_nucleotide_diversity.py \
        --acute /path/to/acute_diversity_summary.csv \
        --chronic /path/to/chronic_diversity_summary.csv \
        --output nucleotide_diversity_boxplot.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import click
from scipy.stats import mannwhitneyu
from matplotlib.ticker import LogFormatterSciNotation
import logging
from pathlib import Path

# === Logging setup ===
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# === Palette and labels ===
PALETTE = {"Acute": "#1f77b4", "Chronic": "#ff7f0e"}
PLOT_TITLE = "Nucleotide Diversity (π) After Replicate Concordance Filtering"
Y_LABEL = "Nucleotide Diversity (π)"
X_LABEL = "Source"

@click.command()
@click.option("--acute", type=click.Path(exists=True, dir_okay=False), required=True, help="CSV file for acute samples")
@click.option("--chronic", type=click.Path(exists=True, dir_okay=False), required=True, help="CSV file for chronic samples")
@click.option("--output", type=click.Path(dir_okay=False), default="nucleotide_diversity_boxplot.png", show_default=True, help="Output PNG filename")
def main(acute, chronic, output):
    # Load data
    logging.info("Loading input data...")
    df_acute = pd.read_csv(acute)
    df_chronic = pd.read_csv(chronic)
    df_acute["source"] = "Acute"
    df_chronic["source"] = "Chronic"
    combined = pd.concat([df_acute, df_chronic], ignore_index=True)

    # Mann–Whitney U test
    logging.info("Performing Mann–Whitney U test...")
    acute_vals = combined[combined["source"] == "Acute"]["pi_merged"].dropna()
    chronic_vals = combined[combined["source"] == "Chronic"]["pi_merged"].dropna()
    u_stat, p_val = mannwhitneyu(acute_vals, chronic_vals, alternative="two-sided")
    logging.info(f"Mann–Whitney U test: U = {u_stat:.2f}, p = {p_val:.2e}")

    # Plotting
    sns.set(style="whitegrid", font_scale=1.4)
    fig, ax = plt.subplots(figsize=(8, 6))

    sns.boxplot(
        data=combined,
        x="source",
        y="pi_merged",
        palette=PALETTE,
        width=0.5,
        linewidth=1.5,
        fliersize=3,
        ax=ax
    )

    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))

    ax.set_xlabel(X_LABEL, fontsize=14)
    ax.set_ylabel(Y_LABEL, fontsize=14)
    ax.set_title(PLOT_TITLE, fontsize=16, pad=20)

    # Determine annotation y position
    box_data = [combined[combined["source"] == grp]["pi_merged"].dropna() for grp in ["Acute", "Chronic"]]
    whisker_tops = [np.percentile(vals, 75) + 1.5 * (np.percentile(vals, 75) - np.percentile(vals, 25)) for vals in box_data]
    outlier_tops = [
        vals[vals > whisker].max() if (vals > whisker).any() else whisker
        for vals, whisker in zip(box_data, whisker_tops)
    ]
    y_star = max(outlier_tops) * 2

    ax.text(0.5, y_star, "***", ha="center", va="bottom", fontsize=18)
    ax.set_ylim(top=y_star * 1.2)

    plt.tight_layout()
    plt.savefig(output, dpi=600, bbox_inches="tight")
    plt.show()
    logging.info(f"Plot saved to: {output}")

if __name__ == "__main__":
    main()
