import logging
import pandas as pd
import re
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple
from mutation_filtering.core.mutation import Mutation
from mutation_filtering.core.filters import apply_soft_filter_df_with_reasons

LAURING_BASE = Path("/sternadi/home/volume3/sars_cov_2/PRJNA889424")
def acute_discover_samples(ivar_dir: Path, output_dir: Path) -> Tuple[List[Tuple], dict]:
    metadata_df = pd.read_csv(LAURING_BASE / "sra_metadata.csv").drop_duplicates()
    metadata_dict = metadata_df.set_index("srr").to_dict("index")
    sample_tasks = []
    
    # Nested defaultdicts: patients[patient][timepoint][rep] = srr_dir
    patients = defaultdict(lambda: defaultdict(dict))

    for srr_dir in ivar_dir.iterdir():
        if not srr_dir.is_dir():
            continue

        srr = srr_dir.name
        if not srr.startswith("SRR"):
            logging.warning(f"Skipping invalid SRR name: {srr}")
            continue

        sample_name = metadata_dict[srr]["sample_name"]
        patient, rep_label = sample_name.split("-")[0], sample_name.split("-")[1]
        rep = "rep1" if rep_label == "A" else "rep2"

        patients[patient]["0"][rep] = srr_dir
        sample_out_dir = output_dir / patient / "0" / rep
        sample_tasks.append((sample_name, srr_dir, sample_out_dir, patient, "0", rep))

    return sample_tasks, {k: dict(v) for k, v in patients.items()}

def chronic_discover_samples(ivar_dir: Path, output_dir: Path) -> Tuple[List[Tuple], dict]:
    sample_tasks = []
    patients = {}

    for patient_dir in ivar_dir.iterdir():
        if not patient_dir.is_dir():
            continue

        patient = patient_dir.name
        if not re.fullmatch(r"[PN]\d+", patient):
            logging.warning(f"Skipping invalid patient name: {patient}")
            continue

        patients[patient] = {}
        for tp_dir in patient_dir.iterdir():
            if not tp_dir.is_dir():
                continue
            timepoint = tp_dir.name
            if not re.fullmatch(r"\d+", timepoint):
                logging.warning(f"Skipping invalid timepoint: {patient}/{timepoint}")
                continue

            patients[patient][timepoint] = {}
            for rep_dir in tp_dir.iterdir():
                if not rep_dir.is_dir():
                    continue
                rep = rep_dir.name
                if rep not in {"rep1", "rep2"}:
                    logging.warning(f"Skipping invalid replicate: {patient}/{timepoint}/{rep}")
                    continue

                sample_name = f"{patient}_{timepoint}_{rep}"
                sample_out_dir = output_dir / patient / timepoint / rep
                sample_tasks.append((sample_name, rep_dir, sample_out_dir, patient, timepoint, rep))

    return sample_tasks, patients

def parse_ivar_file(filepath: Path) -> List[Mutation]:
    try:
        df = pd.read_csv(filepath, sep='\t')
    except pd.errors.EmptyDataError:
        logging.warning(f"{filepath} is empty")
        return []

    if df.empty:
        logging.warning(f"{filepath} is empty after reading")
        return []

    applied = df.apply(Mutation.from_row, axis=1)
    if isinstance(applied, pd.DataFrame):
        raise RuntimeError(f"{filepath} -> apply() returned a DataFrame, not Series. Check Mutation.from_row.")

    return applied.tolist()


def process_sample(
    sample_info: Tuple[str, Path, Path, str, str, str],
    suspicious_set: set,
    problematic_set: set
):
    sample_name, rep_dir, sample_out_dir, patient, timepoint, rep = sample_info

    masked_files = list(rep_dir.rglob("*.masked.tsv"))
    if not masked_files:
        logging.warning(f"No variant files found for {sample_name}")
        return

    mutations = []
    for file in masked_files:
        mutations.extend(parse_ivar_file(file))

    if not mutations:
        logging.warning(f"No mutations for sample {sample_name}")
        return

    sample_out_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{patient}_{timepoint}_{rep}"

    pd.DataFrame([m.to_dict() for m in mutations if m.ALT_FREQ > 0.01]).to_csv(sample_out_dir / f"{prefix}_unfiltered.csv", index=False)

    passed, failed_df = apply_soft_filter_df_with_reasons(mutations, suspicious_set, problematic_set)

    pd.DataFrame([m.to_dict() for m in passed]).to_csv(sample_out_dir / f"{prefix}_soft.csv", index=False)
    failed_df.to_csv(sample_out_dir / f"{prefix}_dropped_soft.csv", index=False)

    logging.info(f"{sample_name}: {len(mutations)} total, {len(passed)} passed soft filter")

    return {
        "sample": sample_name,
        "unfiltered": len(mutations),
        "soft_filtered": len(passed)
    }
