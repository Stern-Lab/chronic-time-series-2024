"""
Codon Position Enrichment Analysis in SARS-CoV-2 Mutations

This script analyzes the distribution of mutations across codon positions (1st, 2nd, 3rd) in coding regions
of SARS-CoV-2, normalized by reference base composition. It processes filtered iVar mutation calls and outputs:
- Normalized mutation proportions by frequency thresholds and bins
- Statistical enrichment tests using the Binomial test
- Publication-ready plots
- CSV of the raw data used in plotting

Usage:
    python analyze_codon_position_enrichment.py \
        --gff path/to/annotations.gff.gz \
        --fasta path/to/reference.fasta \
        --input-dir path/to/ivar_output_directory \
        --output-plot path/to/output_plot.png

Expected input:
    - GFF file with CDS annotations (compressed with gzip)
    - FASTA reference genome
    - iVar output directory containing `masked_variants.tsv` files

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
from collections import Counter
from Bio import SeqIO
from pathlib import Path
from scipy.stats import binomtest
import click

# === Constants ===
palette_codon = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c"}
sns.set(style="whitegrid", font_scale=1.3)

# === Helper Functions ===
def load_reference_positions(fasta_path, positions_of_interest):
    ref_bases = {}
    record = SeqIO.read(fasta_path, "fasta")
    for pos in positions_of_interest:
        ref_bases[pos] = record.seq[pos-1].upper()
    return ref_bases

def load_gene_positions(gff_path):
    entries = []
    with gzip.open(gff_path, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) == 9 and fields[2] == "CDS":
                start, end = int(fields[3]), int(fields[4])
                phase = int(fields[7]) if fields[7].isdigit() else 0
                entries.append((start, end, phase))
    return entries

def assign_codon_positions(gene_entries):
    codon_positions = {}
    for start, end, phase in gene_entries:
        for i, pos in enumerate(range(start, end + 1)):
            codon_positions[pos] = ((i + phase) % 3) + 1
    return codon_positions

def load_mutations(input_dir, codon_positions):
    files = list(input_dir.glob("**/masked_variants.tsv"))
    dfs = []
    for f in files:
        try:
            temp = pd.read_csv(f, sep="\t")
            temp.columns = temp.columns.str.strip()
            dfs.append(temp)
        except Exception:
            pass

    df = pd.concat(dfs, ignore_index=True)
    df = df[
        df["POS"].notnull() &
        df["ALT_FREQ"].notnull() &
        (df["ALT_FREQ"] > 0) & (df["ALT_FREQ"] < 0.5) &
        (df["ALT_QUAL"] >= 30)
    ]
    df["POS"] = df["POS"].astype(int)
    df = df.dropna(subset=["ALT_FREQ"])
    df = df[df["ALT"].apply(lambda x: isinstance(x, str) and len(x) == 1)]
    df["codon_position"] = df["POS"].map(codon_positions)
    return df.dropna(subset=["codon_position"])

def compute_base_counts(ref_bases, codon_positions):
    base_counts = {1: Counter(), 2: Counter(), 3: Counter()}
    for pos, codon_pos in codon_positions.items():
        base = ref_bases.get(pos, "N")
        if base in "ACGT":
            base_counts[codon_pos][base] += 1
    return base_counts

def compute_normalized_rates(sub_df, base_counts):
    mut_counts = {1: Counter(), 2: Counter(), 3: Counter()}
    for _, row in sub_df.iterrows():
        ref_base = row["REF"]
        codon_pos = row["codon_position"]
        if ref_base in "ACGT":
            mut_counts[codon_pos][ref_base] += 1

    norm = {}
    for codon_pos in [1, 2, 3]:
        rate = sum(
            mut_counts[codon_pos][b] / base_counts[codon_pos][b]
            for b in "ACGT" if base_counts[codon_pos][b] > 0
        )
        norm[codon_pos] = rate
    return norm

def normalize_mutations(mut_df, base_counts, thresholds, bin_edges):
    threshold_records = []
    bin_records = []

    for lower, upper in zip(bin_edges[:-1], bin_edges[1:]):
        sub = mut_df[(mut_df["ALT_FREQ"] >= lower) & (mut_df["ALT_FREQ"] < upper)]
        total_mutations = len(sub)
        norm = compute_normalized_rates(sub, base_counts)
        total_norm = sum(norm.values())
        effect_sizes = [(norm[codon_pos] / total_norm) - (1/3) if total_norm > 0 else np.nan for codon_pos in [1, 2, 3]]
        p_values = [binomtest(
            int(norm[codon_pos] * total_mutations / total_norm),
            total_mutations,
            p=1/3,
            alternative="greater"
        ).pvalue if total_mutations > 0 else 1.0 for codon_pos in [1, 2, 3]]

        if total_norm > 0:
            for codon, rate in norm.items():
                bin_records.append({
                    "lower_edge": lower,
                    "codon_position": codon,
                    "normalized_proportion": rate / total_norm,
                    "total_mutations": total_mutations,
                    "p_value": p_values[codon - 1],
                    "effect_size": effect_sizes[codon - 1]
                })

    return pd.DataFrame(threshold_records), pd.DataFrame(bin_records)

def plot_normalized_proportions_bins_only(bin_df, bin_edges, output_plot):
    _, ax = plt.subplots(figsize=(12, 8))

    bin_label_map = {}
    for lower, upper in zip(bin_edges[:-1], bin_edges[1:]):
        subset = bin_df[bin_df["lower_edge"] == lower]
        if not subset.empty:
            n_mut = subset["total_mutations"].iloc[0]
            label = f"[{lower:.1g},{upper:.1g}]\n(n={n_mut})"
            bin_label_map[lower] = label

    bin_df["bin_label"] = bin_df["lower_edge"].map(bin_label_map)
    bin_df["bin_label"] = pd.Categorical(
        bin_df["bin_label"],
        categories=[bin_label_map[x] for x in bin_edges[:-1]],
        ordered=True
    )

    sns.lineplot(
        data=bin_df,
        x="bin_label",
        y="normalized_proportion",
        hue="codon_position",
        hue_order=[1, 2, 3],
        style="codon_position",
        markers=True,
        dashes=False,
        palette=palette_codon,
        ax=ax
    )

    for codon_pos in [1, 2, 3]:
        sub_df = bin_df.query(f"codon_position == {codon_pos}")
        for _, row in sub_df.drop_duplicates("lower_edge").iterrows():
            if (row["p_value"] is not None and row["p_value"] < 0.05 
                and abs(row["effect_size"]) > 0.02):
                ax.text(row["bin_label"], row["normalized_proportion"] + 0.05, "*", 
                        ha='center', va='bottom', color=palette_codon[codon_pos], fontsize=16)

    ax.axhline(1/3, linestyle="--", color="gray", linewidth=1)
    ax.set_ylabel("Normalized Mutation Proportion", fontsize=14, labelpad=10)
    ax.set_xlabel("Mutation Frequency Bin", fontsize=14, labelpad=10)
    ax.set_title(
        "Base Composition-Normalized Mutation Proportions Across Codon Positions by Frequency Bin\n"
        r"(* p < 0.05 and effect size > 2% for enrichment over 1/3; Binomial test)",
        fontsize=16, pad=20
    )

    ax.set_ylim(bin_df["normalized_proportion"].min() * 0.8, bin_df["normalized_proportion"].max() * 1.2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.legend(title="Codon Position", title_fontsize=13, fontsize=12)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_plot, dpi=600)
    print(f"Saved bin-only plot to {output_plot}")

    bin_df_out = bin_df.copy()
    bin_df_out["bin_label"] = bin_df_out["lower_edge"].map(bin_label_map)
    output_csv = output_plot.with_suffix("_raw.csv")
    bin_df_out.to_csv(output_csv, index=False)
    print(f"Saved raw plotting data to {output_csv}")

# === Entry Point ===
@click.command()
@click.option('--gff', type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True, help="Path to compressed GFF file (.gff.gz)")
@click.option('--fasta', type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True, help="Path to reference genome FASTA file")
@click.option('--input-dir', type=click.Path(exists=True, file_okay=False, path_type=Path), required=True, help="Directory containing masked_variants.tsv files")
@click.option('--output-plot', type=click.Path(dir_okay=False, path_type=Path), required=True, help="Path to save the output plot")
def main(gff, fasta, input_dir, output_plot):
    gene_entries = load_gene_positions(gff)
    positions_of_interest = {pos for start, end, _ in gene_entries for pos in range(start, end+1)}
    ref_bases = load_reference_positions(fasta, positions_of_interest)
    codon_positions = assign_codon_positions(gene_entries)
    mut_df = load_mutations(input_dir, codon_positions)
    base_counts = compute_base_counts(ref_bases, codon_positions)

    thresholds = [0, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]
    bin_edges = [0, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.5]

    _, bin_df = normalize_mutations(mut_df, base_counts, thresholds, bin_edges)
    plot_normalized_proportions_bins_only(bin_df, bin_edges, output_plot)

if __name__ == "__main__":
    main()
