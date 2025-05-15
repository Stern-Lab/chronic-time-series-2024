import pandas as pd
import logging
import re
from pathlib import Path
from collections import defaultdict

def track_dropped_mutations(output_dir: Path) -> pd.DataFrame:
    """
    Traverse output directory and collect all dropped mutations
    from soft and replicate filtering. Outputs a summary CSV.
    """
    dropped_files = list(output_dir.rglob("*_dropped_soft.csv")) + list(output_dir.rglob("*_dropped_merged.csv"))
    mutation_drop_locations = defaultdict(list)

    for file in dropped_files:
        try:
            df = pd.read_csv(file)
            if df.empty or "POS" not in df.columns or "mutation" not in df.columns:
                logging.warning(f"Skipping invalid or empty dropped file: {file}")
                continue

            phase = "soft" if "dropped_soft" in file.name else "replicate"

            # Extract sample metadata
            match = re.search(r"([PN]\d+)/(\d+)/(rep[12])?", str(file))
            patient, timepoint, rep = ("unknown", "unknown", "unknown")
            if match:
                patient, timepoint, rep = match.groups()
                rep = rep or "merged"

            metadata = {
                "patient": patient,
                "timepoint": timepoint,
                "rep": rep,
                "phase": phase
            }

            # Create full record per mutation with reason if available
            for _, row in df.iterrows():
                record = {
                    "mutation": row.get("mutation"),
                    "POS": row.get("POS"),
                    "drop_reason": row.get("drop_reason") if "drop_reason" in df.columns else None,
                    **metadata
                }
                mutation_drop_locations[record["mutation"]].append(record)

        except Exception as e:
            logging.error(f"Error processing dropped file {file}: {e}")

    # Flatten all entries
    all_records = [record for entries in mutation_drop_locations.values() for record in entries]
    result_df = pd.DataFrame(all_records)

    if not result_df.empty:
        result_df = result_df.sort_values(by=["POS", "patient", "timepoint"])
        out_file = output_dir / "dropped_mutation_locations.csv"
        result_df.to_csv(out_file, index=False)
        logging.info(f"✅ Dropped mutation summary written to: {out_file}")
        logging.info(f"🔹 Tracked {len(result_df)} mutation drop events")
    else:
        logging.warning("⚠️ No dropped mutations found to track")

    return result_df
