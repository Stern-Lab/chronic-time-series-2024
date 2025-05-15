#!/usr/bin/env python

"""
Acute SARS-CoV-2 Diversity Analysis Pipeline (CLI)

This script processes per-sample mutation CSVs and depth CSVs from iVar runs,
computes nucleotide diversity (π) and overdispersion, and generates:
  - A CSV summary table of diversity metrics
  - A regression plot vs. Ct values

Expected formats:
- mutation_analysis/<patient>/0/<rep>/<patient>_0_<rep>_<filter>.csv
- ivar_runs/<SRR>/<SRR>.depth.csv
- metadata CSV with: srr, sample_name ("<patient>-A/B"), ct_value

Author: Bar Raveh
Lab: Stern Lab, Tel Aviv University
"""

import logging
from pathlib import Path

import click
import pandas as pd

from diversity.plot_diversity import plot_metrics, prepare_regression_plot_df
from diversity.analyze_diversity import (
    process_mutation_file,
    load_depth_file,
    calculate_pi_diversity,
    calculate_overdispersion,
)

# === Logging ===
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def process_acute_mutation_outputs(
    mutation_analysis_path: Path,
    ivar_runs_path: Path,
    metadata_path: Path,
) -> list[dict]:
    metadata_df = pd.read_csv(metadata_path).drop_duplicates()
    metadata_dict = metadata_df.set_index("srr").to_dict("index")

    all_results = []
    processed_patients = set()

    for srr, row in metadata_dict.items():
        sample_name = row["sample_name"]
        try:
            patient, rep_label = sample_name.split("-")
        except ValueError:
            logger.warning(f"Skipping malformed sample name: {sample_name}")
            continue

        rep = "rep1" if rep_label == "A" else "rep2"
        timepoint = "0"

        if patient in processed_patients:
            continue

        rep_paths = {
            "rep1": {
                "unfiltered": mutation_analysis_path / patient / timepoint / "rep1" / f"{patient}_{timepoint}_rep1_unfiltered.csv",
                "soft": mutation_analysis_path / patient / timepoint / "rep1" / f"{patient}_{timepoint}_rep1_soft.csv",
            },
            "rep2": {
                "unfiltered": mutation_analysis_path / patient / timepoint / "rep2" / f"{patient}_{timepoint}_rep2_unfiltered.csv",
                "soft": mutation_analysis_path / patient / timepoint / "rep2" / f"{patient}_{timepoint}_rep2_soft.csv",
            },
        }
        merged_path = mutation_analysis_path / patient / timepoint / f"{patient}_{timepoint}_merged.csv"

        rep_srrs = {
            "rep1": next((s for s, r in metadata_dict.items() if r["sample_name"] == f"{patient}-A"), None),
            "rep2": next((s for s, r in metadata_dict.items() if r["sample_name"] == f"{patient}-B"), None),
        }

        for rep in ["rep1", "rep2"]:
            srr = rep_srrs[rep]
            for category, path in rep_paths[rep].items():
                if not path.exists():
                    continue
                depth_path = ivar_runs_path / srr / f"{srr}.depth.csv" if srr else None
                depth_df = load_depth_file(depth_path) if depth_path and depth_path.exists() else None
                if depth_df is None:
                    continue

                df = process_mutation_file(path, depth_df)
                if df is None:
                    continue

                pi = calculate_pi_diversity(df, f"[{path}]")
                od = calculate_overdispersion(df, f"[{path}]")

                all_results.append({
                    "patient": patient,
                    "timepoint": timepoint,
                    "replicate": rep,
                    "category": category,
                    "sample_name": f"{patient}_{timepoint}_{rep}_{category}",
                    "pi_value": pi,
                    "overdispersion": od,
                })

        # Merged
        if merged_path.exists() and rep_srrs["rep1"] and rep_srrs["rep2"]:
            d1 = load_depth_file(ivar_runs_path / rep_srrs["rep1"] / f"{rep_srrs['rep1']}.depth.csv")
            d2 = load_depth_file(ivar_runs_path / rep_srrs["rep2"] / f"{rep_srrs['rep2']}.depth.csv")
            if d1 is not None and d2 is not None:
                depth_df = pd.concat([d1, d2]).groupby("POS")["depth"].sum().reset_index()
                df = process_mutation_file(merged_path, depth_df, is_merged=True)
                if df is not None:
                    pi = calculate_pi_diversity(df, f"[{merged_path}]")
                    od = calculate_overdispersion(df, f"[{merged_path}]")
                    all_results.append({
                        "patient": patient,
                        "timepoint": timepoint,
                        "replicate": "merged",
                        "category": "merged",
                        "sample_name": f"{patient}_{timepoint}_merged",
                        "pi_value": pi,
                        "overdispersion": od,
                    })

        processed_patients.add(patient)

    return all_results


@click.command()
@click.argument("mutation_analysis_path", type=click.Path(exists=True, path_type=Path))
@click.argument("ivar_runs_path", type=click.Path(exists=True, path_type=Path))
@click.argument("metadata_path", type=click.Path(exists=True, path_type=Path))
@click.argument("output_dir", type=click.Path(file_okay=False, path_type=Path))
@click.option("--project-name", default="acute", help="Project name prefix for output files.")
def analyze_acute(
    mutation_analysis_path: Path,
    ivar_runs_path: Path,
    metadata_path: Path,
    output_dir: Path,
    project_name: str,
):
    """Compute nucleotide diversity and overdispersion for acute SARS-CoV-2 samples."""
    logger.info("Starting acute diversity analysis...")
    results = process_acute_mutation_outputs(mutation_analysis_path, ivar_runs_path, metadata_path)

    if not results:
        logger.error("No results to analyze.")
        return

    df = pd.DataFrame(results)

    df["label"] = df.apply(
        lambda row: f"{row['replicate']}_{row['category']}" if row["replicate"] != "merged" else "merged",
        axis=1
    )
    wide_df = df.pivot_table(
        index=["patient", "timepoint"],
        columns="label",
        values=["pi_value", "overdispersion"]
    )
    wide_df.columns = [f"pi_{c}" if m == "pi_value" else f"od_{c}" for m, c in wide_df.columns]
    wide_df.reset_index(inplace=True)
    wide_df = wide_df.reindex(sorted(wide_df.columns), axis=1)

    output_csv = output_dir / f"{project_name}_diversity_summary.csv"
    wide_df.to_csv(output_csv, index=False)
    logger.info(f"Saved diversity summary to {output_csv}")

    ct_df = pd.read_csv(metadata_path)

    regression_df = prepare_regression_plot_df(
        wide_df,
        ct_df=ct_df,
        ct_column="ct_value",
        ct_sample_column="sample_name",
        chronic=False,
    )

    output_plot_path = output_dir / f"{project_name}_pi_diversity_plot.png"
    plot_metrics(
        regression_df,
        output_path=output_plot_path,
        ct_column="ct_value",
        chronic=False,
    )

    logger.info("Acute diversity analysis complete.")


if __name__ == "__main__":
    analyze_acute()
