import pandas as pd
import logging
from pathlib import Path
from typing import Optional, List


def apply_replicate_concordance(
    output_dir: Path,
    patient: str,
    timepoint: str,
    strict: bool = False,
    skipped: Optional[List[str]] = None
) -> Optional[dict]:
    """
    Merge rep1 and rep2 soft-filtered mutation files and apply replicate concordance filter.
    Saves consensus and dropped results to output_dir.
    """
    rep1_file = output_dir / patient / timepoint / "rep1" / f"{patient}_{timepoint}_rep1_soft.csv"
    rep2_file = output_dir / patient / timepoint / "rep2" / f"{patient}_{timepoint}_rep2_soft.csv"

    if not rep1_file.exists() or not rep2_file.exists():
        msg = f"Missing soft-filtered files for {patient} timepoint {timepoint}"
        _skip_or_raise(msg, f"{patient}_{timepoint}", skipped, strict)
        return

    if rep1_file.stat().st_size == 0 or rep2_file.stat().st_size == 0:
        msg = f"One of the replicate files is empty for {patient} timepoint {timepoint}"
        _skip_or_raise(msg, f"{patient}_{timepoint}", skipped, strict)
        return

    try:
        df1 = pd.read_csv(rep1_file)
        df2 = pd.read_csv(rep2_file)
    except pd.errors.EmptyDataError as e:
        msg = f"Could not parse replicate files for {patient} timepoint {timepoint}: {e}"
        _skip_or_raise(msg, f"{patient}_{timepoint}", skipped, strict)
        return

    merged = pd.merge(df1, df2, on=["POS", "mutation"], how="outer", suffixes=("_x", "_y"))
    merged[["final_freq", "drop_reason"]] = merged.apply(
        lambda row: pd.Series(_replicate_filter_decision_verbose(row)),
        axis=1
    )

    # Additional metrics
    merged["new_freq_x"] = merged["ALT_FREQ_x"].fillna(-1)
    merged["new_freq_y"] = merged["ALT_FREQ_y"].fillna(-1)
    merged["tot_cov"] = merged[["TOTAL_DP_x", "TOTAL_DP_y"]].sum(axis=1)
    merged["tot_base_count"] = merged[["ALT_DP_x", "ALT_DP_y"]].sum(axis=1)

    out_dir = output_dir / patient / timepoint
    out_dir.mkdir(parents=True, exist_ok=True)

    consensus = merged[merged["final_freq"] > 0]
    dropped = merged[merged["final_freq"] <= 0]

    consensus.to_csv(out_dir / f"{patient}_{timepoint}_merged.csv", index=False)
    dropped.to_csv(out_dir / f"{patient}_{timepoint}_dropped_merged.csv", index=False)

    logging.info(
        f"{patient} tp {timepoint}: {len(consensus)} passed replicate filter, {len(dropped)} dropped"
    )

    return {
        "patient": patient,
        "timepoint": timepoint,
        "soft_total": len(df1) + len(df2),
        "replicate_filtered": len(consensus),
    }


def _replicate_filter_decision_verbose(row) -> tuple[float, Optional[str]]:
    f1 = row.get("ALT_FREQ_x")
    f2 = row.get("ALT_FREQ_y")
    d1 = row.get("ALT_DP_x")
    d2 = row.get("ALT_DP_y")

    if pd.isna(f1) or pd.isna(f2) or pd.isna(d1) or pd.isna(d2):
        return -1, "missing_data"
    if f1 == 0 and f2 == 0:
        return 0.0, None
    if (f1 == 0 and f2 > 0) or (f2 == 0 and f1 > 0):
        return -1, "only_one_detected"

    diff = abs(f1 - f2)
    if f1 < 0.5 and f2 < 0.5:
        if diff <= 0.1:
            return (f1 * d1 + f2 * d2) / (d1 + d2), None
        return -1, "low_freq_mismatch"
    if (f1 >= 0.5 or f2 >= 0.5):
        if diff <= 0.3:
            return (f1 * d1 + f2 * d2) / (d1 + d2), None
        return -1, "high_freq_mismatch"

    return -1, "uncategorized"



def _skip_or_raise(msg: str, sample_id: str, skipped: Optional[List[str]], strict: bool):
    logging.warning(msg)
    if skipped is not None:
        skipped.append(sample_id)
    if strict:
        raise FileNotFoundError(msg)
