# SARS-CoV-2 Mutation Filtering & Replicate Concordance Pipeline

This pipeline processes **SARS-CoV-2 variant calls from iVar**, applies filtering to remove unreliable mutations, merges replicate samples, and tracks all dropped mutations across patients and timepoints. It is built for robust, large-scale mutation analysis in longitudinal studies.

---

## Pipeline Summary

### Input:

* iVar `masked_variants.tsv` files, one per replicate

### Filtering Steps:

1. **Soft Filtering**

   * Low depth: `ALT_DP < 50` or `TOTAL_DP < 100`
   * Mutations from a provided suspicious list
   * Sites in a curated VCF of problematic genome positions
2. **Replicate Concordance Filtering**

   * Compares `rep1` and `rep2` at each timepoint
   * Retains mutations with similar support across replicates

### Output:

* Per-replicate mutation tables (unfiltered and soft-filtered)
* Merged mutation tables per timepoint
* Tables of dropped mutations with annotated reasons
* Summary CSVs aggregating filter statistics across samples

---

## Output Directory Structure

```
output_dir/
├── P1/
│   ├── 0/
│   │   ├── rep1/
│   │   │   ├── P1_0_rep1_unfiltered.csv
│   │   │   ├── P1_0_rep1_soft.csv
│   │   │   └── P1_0_rep1_dropped_soft.csv
│   │   ├── rep2/
│   │   │   └── ...
│   │   ├── P1_0_merged.csv
│   │   └── P1_0_dropped_merged.csv
├── dropped_mutation_locations.csv
├── mutation_filtering_summary.csv
└── ...
```

---

## Running the Pipeline

```bash
python -m analyses.mutation_filtering.cli \
    /path/to/metadata.csv \
    /path/to/ivar_dir \
    /path/to/output_dir \
    --workers 8 \
    --suspicious-file /path/to/suspicious_mutations.csv
```

---

## Output Files

| File                             | Description                             |
| -------------------------------- | --------------------------------------- |
| `*_unfiltered.csv`               | Raw iVar mutation calls                 |
| `*_soft.csv`                     | Passed soft filtering                   |
| `*_dropped_soft.csv`             | Dropped after soft filter (with reason) |
| `*_merged.csv`                   | Concordant mutations between replicates |
| `*_dropped_merged.csv`           | Dropped after replicate filtering       |
| `mutation_filtering_summary.csv` | Per-sample filter stats                 |
| `dropped_mutation_locations.csv` | All dropped mutations with drop reason  |

---

## Required Files

| File                             | Description                                             |
| -------------------------------- | ------------------------------------------------------- |
| `sra_metadata.csv`               | Maps sample names to patient and timepoint              |
| `suspicious_mutations.csv`       | List of known artifacts or suspicious mutations         |
| `problematic_sites_sarsCov2.vcf` | Publicly blacklisted sites (auto-downloaded if missing) |

---

## Notes

* Replicates should be located under: `<patient>/<timepoint>/rep1` and `/rep2`
* Sample filenames must follow the format: `P1_0_rep1_soft.csv`, etc.
* Samples with missing or malformed files will be skipped (unless `--strict` is used)
