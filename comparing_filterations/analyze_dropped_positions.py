"""
Analyze dropped mutation positions from a SARS-CoV-2 variant filtering pipeline.

This script compares filtered mutation sites with known problematic and suspicious positions,
and summarizes the reasons for soft and replicate-based mutation filtering.

Usage:
    python analyze_dropped_positions.py \
        --dropped-csv path/to/dropped_mutation_locations.csv \
        --suspicious-csv path/to/suspected_muts_with_soft_filtering_trim_40.csv \
        --output-dir path/to/output_directory

Expected Formats:
- dropped_csv: CSV with columns [mutation, drop_reason, phase] (e.g. "A123C" in mutation).
- suspicious_csv: CSV with a column [mut_nuc] with mutation strings (e.g. "A123C").

Output:
- Venn diagram of mutation positions
- UpSet plot showing overlaps
- Drop reason bar plots for soft and replicate filters
- CSV summaries of drop reasons
"""

import pandas as pd
import matplotlib.pyplot as plt
import re
from pathlib import Path
from venn import venn
from upsetplot import from_indicators, UpSet
import requests
import click

def extract_position(mutation_string: str) -> int:
    match = re.search(r'\d+', str(mutation_string))
    return int(match.group()) if match else None

def load_problematic_positions() -> set:
    vcf_url = "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf"
    vcf_lines = [line for line in requests.get(vcf_url).text.splitlines() if not line.startswith("#")]
    return set(int(line.split("\t")[1]) for line in vcf_lines)

def plot_venn_upset(dropped_df: pd.DataFrame, suspicious_df: pd.DataFrame, output_dir: Path):
    problematic_positions = load_problematic_positions()
    suspicious_df["POS"] = suspicious_df["mut_nuc"].apply(extract_position)
    suspicious_positions = set(suspicious_df["POS"])

    non_prob_susp_df = dropped_df[~dropped_df["drop_reason"].str.contains("problematic|suspicious", na=False)]
    phase_by_pos = non_prob_susp_df.groupby("POS")["phase"].agg(lambda x: set(x)).reset_index()
    phase_by_pos.columns = ["POS", "phase_set"]

    soft_only = set(phase_by_pos[phase_by_pos["phase_set"] == {"soft"}]["POS"])
    soft_or_rep = set(phase_by_pos[phase_by_pos["phase_set"].apply(lambda s: s.issubset({"soft", "replicate"}))]["POS"])

    venn_sets = {
        "Problematic": problematic_positions,
        "Suspicious": suspicious_positions,
        "Filtered (Soft and Replicate)": soft_or_rep
    }
    plt.figure(figsize=(8, 8))
    venn(venn_sets)
    plt.title("Overlap of Mutation Positions")
    plt.tight_layout()
    plt.savefig(output_dir / "mutation_position_overlap_venn3.png")

    all_positions = set().union(problematic_positions, suspicious_positions, soft_only, soft_or_rep)
    upset_df = pd.DataFrame({"POS": sorted(all_positions)})
    upset_df["Problematic"] = upset_df["POS"].isin(problematic_positions)
    upset_df["Suspicious"] = upset_df["POS"].isin(suspicious_positions)
    upset_df["SoftOnly"] = upset_df["POS"].isin(soft_only)
    upset_df["Soft+Replicate"] = upset_df["POS"].isin(soft_or_rep)

    plt.figure(figsize=(10, 6))
    data = from_indicators(upset_df[["Problematic", "Suspicious", "SoftOnly", "Soft+Replicate"]])
    upset = UpSet(data, show_counts=True, subset_size="count")
    upset.plot()
    plt.suptitle("UpSet Plot: Mutation Position Membership")
    plt.tight_layout()
    plt.savefig(output_dir / "mutation_position_upset_4set.png")

def summarize_drop_reasons(df: pd.DataFrame, phase: str, output_dir: Path):
    phase_df = df[df["phase"].str.startswith(phase)].copy()
    reason_series = phase_df["drop_reason"].dropna().str.split(";").explode()
    reason_counts = reason_series.value_counts().sort_values(ascending=False)
    reason_df = reason_counts.reset_index()
    reason_df.columns = ["drop_reason", "count"]
    reason_df.to_csv(output_dir / f"{phase}_drop_reason_counts.csv", index=False)

    plt.figure(figsize=(10, 5))
    plt.bar(reason_df["drop_reason"], reason_df["count"], color="orange" if phase == "soft" else "tomato")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mutation Count")
    plt.title(f"Drop Reasons for {phase.capitalize()} Filtering")
    plt.tight_layout()
    plt.savefig(output_dir / f"{phase}_drop_reason_counts.png")

@click.command()
@click.option('--dropped-csv', type=click.Path(exists=True, path_type=Path), required=True, help='Path to dropped_mutation_locations.csv')
@click.option('--suspicious-csv', type=click.Path(exists=True, path_type=Path), required=True, help='Path to suspicious mutations CSV')
@click.option('--output-dir', type=click.Path(path_type=Path), required=True, help='Directory to save output files')
def main(dropped_csv: Path, suspicious_csv: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    dropped_df = pd.read_csv(dropped_csv)
    dropped_df["POS"] = dropped_df["mutation"].apply(extract_position)
    suspicious_df = pd.read_csv(suspicious_csv)

    plot_venn_upset(dropped_df, suspicious_df, output_dir)
    summarize_drop_reasons(dropped_df, "soft", output_dir)
    summarize_drop_reasons(dropped_df, "replicate", output_dir)

    print("\n✅ Done! Files created in:", output_dir)

if __name__ == "__main__":
    main()