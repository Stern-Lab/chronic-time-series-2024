#!/usr/bin/env python

"""
Analyze SARS-CoV-2 Diversity

This script computes nucleotide Pi diversity (π) and overdispersion across
replicates and merged mutation calls. It summarizes metrics across patients
and timepoints, and correlates them with Ct values.

Inputs:
--mutation-analysis-path: directory with patient/timepoint/replicate CSVs:
    mutation_analysis/
    └── <patient>/<timepoint>/
        ├── rep1/<patient>_<timepoint>_rep1_[unfiltered|soft].csv
        ├── rep2/<patient>_<timepoint>_rep2_[unfiltered|soft].csv
        └── <patient>_<timepoint>_merged.csv

--ivar-runs-path: contains per-replicate depth TSVs:
    ivar_runs/
    └── <patient>/<timepoint>/<rep>/masked_variants_depth.tsv

--ct-values-file: CSV containing 'sample_name' and Ct column (e.g., 'ct_value')

Outputs:
- <project_name>_diversity_summary.csv: pivoted summary of π and overdispersion
- <project_name>_pi_diversity_plot.png: regression of π vs. Ct

Example usage:
python analyze_diversity.py \
    --mutation-analysis-path ./mutation_analysis \
    --ivar-runs-path ./ivar_runs \
    --ct-values-file ./metadata.csv \
    --output-path ./results \
    --project-name chronic_samples
"""

import logging
import os
from typing import Optional

import click
import pandas as pd
import numpy as np
import concurrent.futures

from diversity.plot_diversity import plot_metrics, prepare_regression_plot_df
from diversity.format import FrequencyTableRow
from diversity.format import (
    COL_BASE_COUNT,
    COL_COVERAGE,
    COL_FREQUENCY,
    COL_BASE_RANK,
    COL_READ_BASE,
    COL_REF_POS
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def calculate_overdispersion(data: pd.DataFrame, log_prefix: str = "") -> float:
    if data.empty:
        logger.warning(f"{log_prefix}Empty input data received for overdispersion analysis")
        return np.nan
    data = data[data[COL_FREQUENCY] > 0]
    mean_f = data[COL_FREQUENCY].mean()
    var_f = data[COL_FREQUENCY].var()

    if var_f >= mean_f * (1 - mean_f):
        logger.warning(f"{log_prefix}Variance too high to estimate Beta-Binomial parameters")
        return np.nan

    α = mean_f * ((mean_f * (1 - mean_f)) / var_f - 1)
    β = (1 - mean_f) * ((mean_f * (1 - mean_f)) / var_f - 1)
    if α <= 0 or β <= 0:
        logger.warning(f"{log_prefix}Invalid α or β estimates")
        return np.nan

    return α / (α + β)


def calculate_pi_diversity(data: pd.DataFrame, log_prefix: str = "") -> float:
    data = data[[COL_REF_POS, COL_READ_BASE, COL_BASE_COUNT, COL_BASE_RANK]].copy()
    data = data[data[COL_BASE_RANK].isin([0, 1])]
    data[COL_BASE_RANK] = data[COL_BASE_RANK].map({0: "Major", 1: "Minor"})

    grouped = data.groupby([COL_REF_POS, COL_BASE_RANK])[COL_BASE_COUNT].agg("max").unstack(COL_BASE_RANK).reset_index()
    if "Minor" not in grouped.columns:
        logger.warning(f"{log_prefix}No minor variants found")
        return 0.0

    grouped["Minor"] = grouped["Minor"].fillna(0)
    grouped["Total"] = grouped["Major"] + grouped["Minor"]
    grouped["pdp"] = grouped.apply(
        lambda row: 0 if row["Minor"] == 0 else (
            row["Total"] * (row["Total"] - 1) -
            ((row["Major"] * (row["Major"] - 1)) + (row["Minor"] * (row["Minor"] - 1)))
        ) / (row["Total"] * (row["Total"] - 1)),
        axis=1
    )
    return grouped["pdp"].mean()


def load_depth_file(depth_path: str) -> Optional[pd.DataFrame]:
    if not os.path.exists(depth_path):
        logger.warning(f"Depth file not found: {depth_path}")
        return None
    try:
        df = pd.read_csv(depth_path, sep="\t", header=None)
        if df.shape[1] == 3:
            df.columns = ["chrom", "POS", "depth"]
        elif df.shape[1] == 2:
            df.columns = ["POS", "depth"]
        else:
            logger.error(f"Unexpected depth file format in {depth_path}. Shape: {df.shape}")
            return None
        return df[["POS", "depth"]]
    except Exception as e:
        logger.error(f"Failed to read depth file {depth_path}: {e}")
        return None


def process_mutation_file(file_path: str, depth_df: pd.DataFrame, is_merged: bool = False) -> Optional[pd.DataFrame]:
    try:
        df = pd.read_csv(file_path)
        mut_df = _extract_merged_mutation_df(df) if is_merged else _extract_single_mutation_df(df)
        if mut_df is None:
            return None

        ref_only_df = _extract_ref_only_df(depth_df, mut_df[COL_REF_POS].unique())
        full_df = pd.concat([mut_df, ref_only_df], ignore_index=True)

        full_df[COL_BASE_RANK] = (
            full_df.groupby(COL_REF_POS)[COL_BASE_COUNT]
            .rank(method="first", ascending=False)
            .sub(1)
        ).astype(pd.Int64Dtype())

        FrequencyTableRow.verify_columns(full_df)
        return full_df

    except Exception as e:
        logger.error(f"Error processing {file_path}: {e}")
        return None


def _extract_single_mutation_df(df: pd.DataFrame) -> Optional[pd.DataFrame]:
    required = ["POS", "REF", "ALT", "ALT_DP", "TOTAL_DP", "ALT_FREQ"]
    if missing := [c for c in required if c not in df.columns]:
        logger.error(f"Missing columns in single-replicate file: {missing}")
        return None

    df = df[df["ALT"].isin(["A", "C", "G", "T"])]
    alt_df = pd.DataFrame({
        COL_REF_POS: df["POS"],
        COL_READ_BASE: df["ALT"],
        COL_BASE_COUNT: df["ALT_DP"],
        COL_COVERAGE: df["TOTAL_DP"],
        COL_FREQUENCY: df["ALT_FREQ"]
    })
    ref_df = pd.DataFrame({
        COL_REF_POS: df["POS"],
        COL_READ_BASE: df["REF"],
        COL_BASE_COUNT: df["TOTAL_DP"] - df["ALT_DP"],
        COL_COVERAGE: df["TOTAL_DP"],
        COL_FREQUENCY: 1.0
    })
    return pd.concat([alt_df, ref_df], ignore_index=True)


def _extract_merged_mutation_df(df: pd.DataFrame) -> Optional[pd.DataFrame]:
    required = [
        "POS", "ALT_x", "new_freq_x", "ALT_y", "new_freq_y", "final_freq",
        "tot_cov", "tot_base_count", "REF_x", "REF_y"
    ]
    if missing := [c for c in required if c not in df.columns]:
        logger.error(f"Missing columns in merged file: {missing}")
        return None

    df = df[df["ALT_x"].isin(["A", "C", "G", "T"]) | df["ALT_y"].isin(["A", "C", "G", "T"])]
    df["merged_alt"] = np.where(df["new_freq_x"] >= df["new_freq_y"], df["ALT_x"], df["ALT_y"])

    alt_df = pd.DataFrame({
        COL_REF_POS: df["POS"],
        COL_READ_BASE: df["merged_alt"],
        COL_BASE_COUNT: df["tot_base_count"],
        COL_COVERAGE: df["tot_cov"],
        COL_FREQUENCY: df["final_freq"]
    })
    ref_base = df["REF_x"].fillna(df["REF_y"])
    ref_df = pd.DataFrame({
        COL_REF_POS: df["POS"],
        COL_READ_BASE: ref_base,
        COL_BASE_COUNT: df["tot_cov"] - df["tot_base_count"],
        COL_COVERAGE: df["tot_cov"],
        COL_FREQUENCY: 1.0
    })
    return pd.concat([alt_df, ref_df], ignore_index=True)


def _extract_ref_only_df(depth_df: pd.DataFrame, mutated_pos: np.ndarray) -> pd.DataFrame:
    ref_only = depth_df[~depth_df["POS"].isin(mutated_pos)].copy()
    return pd.DataFrame({
        COL_REF_POS: ref_only["POS"],
        COL_READ_BASE: "N",
        COL_BASE_COUNT: ref_only["depth"],
        COL_COVERAGE: ref_only["depth"],
        COL_FREQUENCY: 1.0,
        COL_BASE_RANK: 0
    })


def process_one_timepoint(mutation_analysis_path: str, ivar_runs_path: str, patient: str, tp: str) -> list[dict]:
    results = []
    tp_path = os.path.join(mutation_analysis_path, patient, tp)
    rep_dirs = {
        "rep1": os.path.join(tp_path, "rep1"),
        "rep2": os.path.join(tp_path, "rep2"),
    }

    for rep, rep_dir in rep_dirs.items():
        for category in ["unfiltered", "soft"]:
            stem = f"{patient}_{tp}_{rep}_{category}"
            file_path = os.path.join(rep_dir, f"{stem}.csv")
            depth_path = os.path.join(ivar_runs_path, patient, tp, rep, "masked_variants_depth.tsv")
            if not os.path.exists(file_path):
                continue
            depth_df = load_depth_file(depth_path)
            if depth_df is None:
                continue
            df = process_mutation_file(file_path, depth_df)
            if df is None:
                continue
            pi = calculate_pi_diversity(df, f"[{file_path}]")
            od = calculate_overdispersion(df, f"[{file_path}]")
            results.append({
                "patient": patient,
                "timepoint": tp,
                "replicate": rep,
                "category": category,
                "sample_name": stem,
                "pi_value": pi,
                "overdispersion": od,
            })

    # Merged
    merged_file = os.path.join(tp_path, f"{patient}_{tp}_merged.csv")
    rep1_depth = os.path.join(ivar_runs_path, patient, tp, "rep1", "masked_variants_depth.tsv")
    rep2_depth = os.path.join(ivar_runs_path, patient, tp, "rep2", "masked_variants_depth.tsv")

    if os.path.exists(merged_file):
        d1 = load_depth_file(rep1_depth)
        d2 = load_depth_file(rep2_depth)
        if d1 is not None and d2 is not None:
            depth_df = pd.concat([d1, d2]).groupby("POS")["depth"].sum().reset_index()
            df = process_mutation_file(merged_file, depth_df, is_merged=True)
            if df is not None:
                pi = calculate_pi_diversity(df, f"[{merged_file}]")
                od = calculate_overdispersion(df, f"[{merged_file}]")
                results.append({
                    "patient": patient,
                    "timepoint": tp,
                    "replicate": "merged",
                    "category": "merged",
                    "sample_name": f"{patient}_{tp}_merged",
                    "pi_value": pi,
                    "overdispersion": od,
                })

    return results


def run_one_timepoint(args): return process_one_timepoint(*args)


@click.command()
@click.option('--mutation-analysis-path', required=True, type=click.Path(exists=True))
@click.option('--ivar-runs-path', required=True, type=click.Path(exists=True))
@click.option('--output-path', required=True, type=click.Path(exists=True))
@click.option('--ct-values-file', required=True, type=click.Path(exists=True))
@click.option('--ct-column-name', default="ct_value")
@click.option('--ct-sample-column', default="sample_name")
@click.option('--project-name', default=None)
def main(mutation_analysis_path, ivar_runs_path, output_path, ct_values_file,
         ct_column_name, ct_sample_column, project_name):
    allowed_timepoints = {
        "N1": ["0", "45", "165", "189"], "N2": ["0", "91", "105"],
        "N3": ["0", "5",  "30",  "76"], "N4": ["0", "27"],
        "N7": ["0", "21"], "N8": ["0", "13", "16", "19", "23", "26", "31", "37"],
        "P3": ["0", "16", "43", "52", "80"], "P4": ["0", "10", "18", "25", "28"],
        "P5": ["0", "12", "20", "30", "40", "68"],
    }

    pairs = []
    for patient in os.listdir(mutation_analysis_path):
        if patient not in allowed_timepoints:
            continue
        for tp in os.listdir(os.path.join(mutation_analysis_path, patient)):
            if tp in allowed_timepoints[patient]:
                pairs.append((mutation_analysis_path, ivar_runs_path, patient, tp))

    logger.info(f"Processing {len(pairs)} patient/timepoints in parallel")
    all_results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in executor.map(run_one_timepoint, pairs):
            all_results.extend(result)

    if not all_results:
        logger.warning("No valid results found.")
        return

    df = pd.DataFrame(all_results)
    df["label"] = df.apply(lambda row: f"{row['replicate']}_{row['category']}" if row["replicate"] != "merged" else "merged", axis=1)
    wide_df = df.pivot_table(index=["patient", "timepoint"], columns="label", values=["pi_value", "overdispersion"])
    wide_df.columns = [f"pi_{c}" if m == "pi_value" else f"od_{c}" for m, c in wide_df.columns]
    wide_df.reset_index(inplace=True)
    wide_df = wide_df.reindex(sorted(wide_df.columns), axis=1)

    if project_name is None:
        project_name = os.path.basename(mutation_analysis_path.rstrip("/"))
    output_csv = os.path.join(output_path, f"{project_name}_diversity_summary.csv")
    wide_df.to_csv(output_csv, index=False)
    logger.info(f"Saved summary to {output_csv}")

    ct_df = pd.read_csv(ct_values_file)
    regression_df = prepare_regression_plot_df(
        wide_df,
        ct_df=ct_df,
        ct_column=ct_column_name,
        ct_sample_column=ct_sample_column,
    )
    plot_metrics(
        regression_df,
        output_path=os.path.join(output_path, f"{project_name}_pi_diversity_plot.png"),
        ct_column=ct_column_name,
    )
    logger.info("Analysis complete.")


if __name__ == "__main__":
    main()
