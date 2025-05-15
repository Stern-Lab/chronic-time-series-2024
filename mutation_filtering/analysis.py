#!/usr/bin/env python

"""
Batch mutation analysis pipeline for acute SARS-CoV-2 samples.

This script processes iVar mutation outputs, applies filtering strategies,
performs replicate concordance filtering, and summarizes dropped mutations.

Author: Bar Raveh
Lab: Adi Stern Lab, Tel Aviv University
"""

import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List

from mutation_filtering.core.io_utils import (
    acute_discover_samples,
    process_sample
)
from mutation_filtering.core.replicate_filter import apply_replicate_concordance
from mutation_filtering.core.drop_tracker import track_dropped_mutations
from mutation_filtering.core.filters import (
    load_problematic_sites,
    load_suspicious_mutations
)


def run_analysis(
    ivar_dir: str,
    output_dir: str,
    workers: int = 4,
    log_level: str = "INFO",
    strict: bool = False,
    suspicious_file: str = "suspicious_mutations.csv"
) -> None:
    """
    Run the full acute mutation analysis pipeline.

    Args:
        ivar_dir (str): Path to directory containing iVar sample outputs.
        output_dir (str): Path to output directory for filtered results.
        workers (int): Number of threads to use for processing samples.
        log_level (str): Logging level (e.g., 'INFO', 'DEBUG').
        strict (bool): Whether to apply strict replicate concordance filtering.
        suspicious_file (str): Path to TSV or CSV file with suspicious mutations.
    """
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        level=getattr(logging, log_level.upper())
    )

    ivar_path = Path(ivar_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load known problematic and suspicious mutations
    suspicious_df = load_suspicious_mutations(Path(suspicious_file))
    suspicious_set = set(suspicious_df["mutation"])
    problematic_set = load_problematic_sites()

    # Discover all acute samples
    sample_tasks, patients = acute_discover_samples(ivar_path, output_path)
    logging.info(f"Discovered {len(sample_tasks)} samples across {len(patients)} patients.")

    sample_summaries: List[dict] = []

    # === Process each sample in parallel ===
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                process_sample,
                task,
                suspicious_set,
                problematic_set
            ): task[0]
            for task in sample_tasks
        }

        for future in as_completed(futures):
            sample_name = futures[future]
            try:
                result = future.result()
                if result:
                    sample_summaries.append(result)
            except Exception as e:
                logging.error(f"❌ Failed to process sample {sample_name}: {e}")

    # === Apply replicate concordance filtering ===
    skipped = []
    replicate_summary = []

    for patient in patients:
        for timepoint in patients[patient]:
            result = apply_replicate_concordance(
                output_path,
                patient,
                timepoint,
                strict=strict,
                skipped=skipped
            )
            if result:
                replicate_summary.append(result)

    if skipped:
        logging.warning(f"⚠️ Skipped {len(skipped)} timepoints due to missing or invalid replicates: {skipped}")

    # === Track dropped mutations ===
    track_dropped_mutations(output_path)

    logging.info("✅ Mutation analysis completed successfully.")
