"""
Plot SARS-CoV-2 Diversity Metrics

This module provides plotting functions for visualizing nucleotide diversity (π),
overdispersion, and their relationship to clinical metadata (e.g., Ct values) over time.

Input:
- summary_df: Pivoted output from diversity analysis with columns:
    ["patient", "timepoint", "pi_<replicate>_<filter>", "od_<replicate>_<filter>", ...]

- ct_df: (Optional) Metadata table with:
    - sample_name (e.g., "N1_0_rep1_unfiltered")
    - ct_value

Output:
- Regression plots of π vs. Ct
- Diversity timecourse plots by patient
- Aggregate π trends over time
- Ct trends over time

Usage:
This module is designed to be imported in scripts such as `analyze_diversity.py`.
"""

import os
import logging
import math
from typing import Optional
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import LogFormatterSciNotation
from scipy.stats import t, linregress
import numpy as np
import seaborn as sns
import re


plt.rcParams.update({
    "axes.titlesize": 18,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "legend.title_fontsize": 13,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})

logger = logging.getLogger(__name__)

COLOR_MAP = {
        "unfiltered": "#1f77b4",              # blue
        "soft": "#ff7f0e",                    # orange
        "replicate_concordance": "#2ca02c"   # green
    }

def plot_aggregate_diversity_over_time(df: pd.DataFrame, output_path: str):
    df = df[df["pi"].notnull() & df["timepoint"].notnull()].copy()
    df["timepoint"] = df["timepoint"].astype(float)
    df["pi"] = df["pi"].clip(lower=1e-7)

    plt.figure(figsize=(10, 6))
    sns.set(style="ticks", context="notebook")
    
    for divtype in df["diversity_type"].unique():
        if divtype == 'soft':
            continue
        sub = df[df["diversity_type"] == divtype].copy()
        color = COLOR_MAP.get(divtype, "gray")

        # Compute regression stats
        x = sub["timepoint"].values
        y = np.log10(sub["pi"].values)
        slope, intercept, r_value, p_value, _ = linregress(x, y)
        label = f"{divtype} (b={slope:.1g}, R²={r_value**2:.2f}, p={p_value:.1g})"

        # Scatter plot with linear regression line
        sns.regplot(
            x="timepoint",
            y="pi",
            data=sub,
            scatter=True,
            color=color,
            label=label,
            line_kws={"linewidth": 2},
            ci=95,
            truncate=False
        )

    plt.yscale("log")
    plt.gca().yaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))
    plt.xlabel("Timepoint (days since first positive report)")
    plt.ylabel("Nucleotide Diversity (π, log scale)")
    plt.legend(title="Filtering Method", title_fontsize=11, loc="lower right")
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    plt.close()


def plot_diversity_timecourses_by_patient(
    df: pd.DataFrame,
    output_path: str,
    y_col: str = "pi",
    title: str = "Nucleotide Diversity Trajectories by Patient"
):
    df = df[df[y_col].notnull() & df["timepoint"].notnull()].copy()
    df["timepoint"] = df["timepoint"].astype(float)
    patients = sorted(df["patient"].unique())

    n_pat = len(patients)
    cols = 3
    rows = math.ceil(n_pat / cols)

    fig, axes = plt.subplots(
        rows, cols,
        figsize=(cols * 6, rows * 4),
        sharex=False,
        sharey=True
    )
    axes = axes.flatten()

    for idx, (ax, patient_id) in enumerate(zip(axes, patients)):
        sub = df[df["patient"] == patient_id]

        sns.lineplot(
            data=sub,
            x="timepoint",
            y=y_col,
            hue="diversity_type",
            ax=ax,
            marker="o",
            linewidth=1.5,
            err_style="band",
            errorbar=('ci', 95),
            palette="colorblind"
        )

        # Axis titles and labels
        ax.set_title(f"Patient {patient_id}", loc="center")
        ax.set_xlabel("Timepoint (days since first positive report)")
        ax.set_ylabel("Nucleotide Diversity (π)" if idx % cols == 0 else "")

        # X ticks
        timepoints = sorted(sub["timepoint"].unique())
        ax.set_xticks(timepoints)
        ax.set_xticklabels([str(int(tp)) for tp in timepoints], rotation=30, ha='right')

        # Remove legends per axis
        ax.get_legend().remove()

    # Remove empty subplots
    for ax in axes[n_pat:]:
        fig.delaxes(ax)

    # Shared legend under title
    handles, labels = axes[0].get_legend_handles_labels()
    # Set up title first
        # Set up title first with more space above the subplots
    fig.suptitle(
        "Nucleotide Diversity Trajectories by Patient",
        fontsize=24,
    )

    plt.tight_layout(pad=1.5, rect=[0, 0, 1, 0.97])  # almost full height for subplots


    fig.legend(
        handles, labels,
        title="Filtering Method",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.06),
        ncol=3,
        fontsize=16,
        title_fontsize=18,
        frameon=False
    )
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(output_path, bbox_inches="tight", dpi=600)
    plt.close(fig)


def prepare_regression_plot_df(
    summary_df: pd.DataFrame,
    ct_df: Optional[pd.DataFrame] = None,
    ct_column: str = "ct_value",
    ct_sample_column: str = "sample_name",
    chronic: bool = True,
) -> pd.DataFrame:
    required_pi_cols = [
        "pi_rep1_unfiltered", "pi_rep2_unfiltered",
        "pi_rep1_soft",       "pi_rep2_soft",
        "pi_merged"
    ]
    required_od_cols = [
        "od_rep1_unfiltered", "od_rep2_unfiltered",
        "od_rep1_soft",       "od_rep2_soft",
        "od_merged"
    ]
    required_columns = ["patient", "timepoint"] + required_pi_cols + required_od_cols
    missing = [col for col in required_columns if col not in summary_df.columns]
    if missing:
        raise ValueError(f"Missing expected columns in summary_df: {missing}")

    records = []
    for _, row in summary_df.iterrows():
        base = {
            "patient":   row["patient"],
            "timepoint": row["timepoint"],
        }

        for prefix in ["rep1_unfiltered", "rep2_unfiltered", "rep1_soft", "rep2_soft"]:
            pi_col = f"pi_{prefix}"
            od_col = f"od_{prefix}"
            records.append({
                **base,
                "sample_name":    f"{row['patient']}_{row['timepoint']}_{prefix}",
                "diversity_type": prefix.split("_")[1],
                "pi":             row[pi_col],
                "overdispersion": row[od_col],
            })

        records.append({
            **base,
            "sample_name":    f"{row['patient']}_{row['timepoint']}_merged",
            "diversity_type": "replicate_concordance",
            "pi":             row["pi_merged"],
            "overdispersion": row["od_merged"],
        })

    df = pd.DataFrame(records)
    df["sample_name"] = df["sample_name"].str.strip()

    if ct_df is not None:
        if ct_sample_column not in ct_df.columns or ct_column not in ct_df.columns:
            raise ValueError(
                f"CT file is missing one of the required columns: "
                f"{ct_sample_column}, {ct_column}"
            )

        ct_df = ct_df.copy()

        if chronic:
            df["sample_name"] = df["sample_name"].str.strip().str.lower()
            ct_df[ct_sample_column] = (
                ct_df[ct_sample_column]
                    .astype(str)
                    .str.lower()
                    .str.strip()
            )

            # for merged samples, borrow the rep1 CT key
            # e.g. "n1_0_merged" -> "n1_0_rep1"
            df["ct_merge_key"] = (
                df["sample_name"]
                .str.replace(r"_merged$", "_rep1", regex=True)
                .str.extract(r"^(.+?_rep\d)", expand=False)
            )

            # prepare ct lookup
            ct_lookup = (
                ct_df[[ct_sample_column, ct_column]]
                .rename(columns={ct_sample_column: "ct_merge_key"})
            )

        else:
            def sample_to_ct_key(name: str) -> str:
                name = name.strip().lower()
                if name.endswith("_merged"):
                    name = name.replace("_merged", "_rep1")
                match = re.match(r"^(?P<patient>.+?)_.*_(?P<rep>rep\d)(?:_.*)?$", name)
                if not match:
                    return ""
                rep = match.group("rep")
                rep_code = {"rep1": "A", "rep2": "B"}.get(rep, "")
                return f"{match.group('patient')}-{rep_code}".lower()

            df["ct_merge_key"] = df["sample_name"].apply(sample_to_ct_key)
            ct_df["ct_merge_key"] = ct_df[ct_sample_column].astype(str).str.strip().str.lower()
            ct_lookup = ct_df[["ct_merge_key", ct_column]]

        df = df.merge(ct_lookup, on="ct_merge_key", how="left")
        merged_count = df[ct_column].notnull().sum()
        df.drop(columns=["ct_merge_key"], inplace=True)

        logger.info(f"Merged CT values: {merged_count} out of {len(df)}")
        if merged_count == 0:
            logger.warning(
                "⚠️ No CT values matched after merge. "
                "Please check sample naming convention."
            )

    return df


def regression_panel(df, ax, x, y, title, log_scale: bool = True):
    sub = df[[x, y, "diversity_type"]].dropna()
    sub = sub[sub[y] > 0]

    if sub.empty:
        ax.set_title(f"No data for {title}")
        return

    handles = []
    labels = []

    for divtype in ["unfiltered", "soft", "replicate_concordance"]:
        subset = sub[sub["diversity_type"] == divtype]
        if len(subset) < 3:
            continue

        xdata = subset[x].values
        ydata = subset[y].values

        if log_scale:
            ydata = np.clip(ydata, 1e-7, None)
            y_fit_data = np.log10(ydata)
        else:
            y_fit_data = ydata

        # Linear regression
        slope, intercept, r, p, _ = linregress(xdata, y_fit_data)

        x_fit = np.linspace(xdata.min(), xdata.max(), 200)
        y_fit = intercept + slope * x_fit

        # Confidence interval
        n = len(xdata)
        x_mean = np.mean(xdata)
        SE_y = np.sqrt(np.sum((y_fit_data - (intercept + slope * xdata)) ** 2) / (n - 2))
        t_val = t.ppf(0.975, df=n - 2)

        conf_margin = t_val * SE_y * np.sqrt(
            1 / n + ((x_fit - x_mean) ** 2 / np.sum((xdata - x_mean) ** 2))
        )

        y_lower = y_fit - conf_margin
        y_upper = y_fit + conf_margin

        if log_scale:
            y_fit = 10 ** y_fit
            y_lower = 10 ** y_lower
            y_upper = 10 ** y_upper

        # Plot regression line and CI band
        curve = ax.plot(x_fit, y_fit, color=COLOR_MAP[divtype])[0]
        ax.fill_between(x_fit, y_lower, y_upper, alpha=0.2, color=COLOR_MAP[divtype])

        # Scatter
        sns.scatterplot(
            data=subset,
            x=x,
            y=y,
            ax=ax,
            s=30,
            alpha=0.6,
            color=COLOR_MAP[divtype]
        )

        label = f"{divtype} (R²={r**2:.2f}, p={p:.1g})"
        handles.append(curve)
        labels.append(label)

    # Axes labels and formatting
    ax.set_xlabel("Ct Value")
    ax.set_ylabel("Nucleotide Diversity (π, log scale)" if log_scale else "Nucleotide Diversity (π)")

    if log_scale:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(LogFormatterSciNotation())

    ax.legend(
        handles, labels,
        title="Filtering Method",
        fontsize=9,
        title_fontsize=10,
        loc="upper right"
    )
    ax.set_title(title)

def plot_ct_over_time(ct_df: pd.DataFrame, output_path: str,
                      ct_column: str = "ct_value",
                      time_column: str = "timepoint",
                      patient_column: str = "patient",
                      max_time: Optional[int] = None):
    # === Input Validation ===
    required_cols = [ct_column, time_column, patient_column]
    for col in required_cols:
        if col not in ct_df.columns:
            raise ValueError(f"Missing column: {col}")

    df = ct_df.copy()
    df[time_column] = pd.to_numeric(df[time_column], errors="coerce")
    
    df[ct_column] = pd.to_numeric(df[ct_column], errors="coerce")
    df = df.dropna(subset=[ct_column, time_column])

    if max_time:
        df = df[df[time_column] < max_time]

    # === Compute Regression ===
    reg = linregress(df[time_column], df[ct_column])
    slope, intercept, r_value, p_value = reg.slope, reg.intercept, reg.rvalue, reg.pvalue
    annotation = f"b = {slope:.2f}\nR² = {r_value**2:.2f}\np = {p_value:.2g}"

    # === Plot ===
    plt.figure(figsize=(10, 6))
    sns.set(style="whitegrid", context="notebook")

    # Scatterplot
    sns.scatterplot(
        data=df,
        x=time_column,
        y=ct_column,
        hue=patient_column,
        palette="tab10",
        s=60,
        edgecolor="black"
    )

    # Regression line
    sns.regplot(
        data=df,
        x=time_column,
        y=ct_column,
        scatter=False,
        color="black",
        line_kws={"linewidth": 3},
        ci=90
    )

    # Invert y-axis (lower Ct = higher viral load)
    plt.gca().invert_yaxis()

    # Text box with regression stats (top-left)
    plt.gca().text(
        0.02, 0.98, annotation,
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment='top',
        horizontalalignment='left',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", alpha=0.9)
    )

    # Labels and title
    plt.title("Ct Value Over Time Across Patients", fontsize=14)
    plt.xlabel("Timepoint (days since first positive report)", fontsize=13)
    plt.ylabel("Ct Value", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Tighter y-axis if helpful
    plt.ylim(35, 15)  # inverted axis

    # Save
    plt.tight_layout(pad=1.5)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()



def plot_metrics(
    df: pd.DataFrame,
    output_path: str,
    ct_column: str,
    chronic=True,
):
    df = df.copy()
    df["pi"] = pd.to_numeric(df["pi"], errors="coerce")
    df["overdispersion"] = pd.to_numeric(df["overdispersion"], errors="coerce")

    if ct_column not in df.columns:
        logger.error(f"CT column '{ct_column}' not found in dataframe columns: {df.columns.tolist()}")
        return

    _, ax = plt.subplots(1, 1, figsize=(6.7, 5.5))
    regression_panel(df, ax, ct_column, "pi", "Pi Diversity vs CT")

    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    logger.info(f"Saved regression plots to {output_path}")
    plt.close()
    if chronic:
        # ---- Combined per‑patient grid of timecourses ---- #
        combined = df[df["pi"].notnull() & df["timepoint"].notnull()].copy()
        combined["timepoint"] = combined["timepoint"].astype(float)

        corrections = {
            "N1" : 18,
            "N2": 29,
            "N3": 60,
            "N4": 0,
            "N7": 41,
            "N8": 4,
            "P3": 30,
            "P4": 26,
            "P5": 5, 
        }
        combined["timepoint"] = combined.apply(
            lambda row: row["timepoint"] + corrections.get(row["patient"], 0),
            axis=1
        )

        print(combined.head(10))

        diversity_timecourses_by_patient_output = os.path.splitext(output_path)[0] + "_per_patient_grid_timecourse.png"
        plot_diversity_timecourses_by_patient(
            df=combined,
            output_path=diversity_timecourses_by_patient_output
        )

        logger.info(f"Saved diversity-over-time plot: {diversity_timecourses_by_patient_output}")

        ct_plot_path = os.path.splitext(output_path)[0] + "_ct_over_time.png"
        try:
            plot_ct_over_time(combined, output_path=ct_plot_path)
            logger.info(f"Saved Ct-over-time plot: {ct_plot_path}")
        except Exception as e:
            logger.warning(f"Could not generate Ct-over-time plot: {e}")
        
        plot_aggregate_diversity_over_time(combined, os.path.splitext(output_path)[0] + "_aggregated_timecourse.png")




