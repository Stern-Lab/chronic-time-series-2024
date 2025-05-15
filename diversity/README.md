### `diversity/` – Nucleotide Diversity & Overdispersion Analysis

This module computes nucleotide diversity (π) and overdispersion metrics from mutation call tables (e.g., iVar outputs), and generates summary plots for viral evolution analysis.

---

## Contents

| Script                 | Description                                                                                                                         |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `analyze_diversity.py` | Computes π and overdispersion for each sample (unfiltered, soft-filtered, merged). Outputs a summary CSV and diversity vs. Ct plot. |
| `plot_diversity.py`    | Generates plots of π over time and against Ct values, including per-patient timecourses and regression analyses.                    |
| `format.py`            | Internal column name conventions and validation logic for mutation tables.                                                          |

---

## Expected Input Structure

### 1. Mutation Analysis Directory

```
mutation_analysis/
└── <patient>/<timepoint>/
    ├── rep1/<patient>_<timepoint>_rep1_unfiltered.csv
    ├── rep1/<patient>_<timepoint>_rep1_soft.csv
    ├── rep2/<patient>_<timepoint>_rep2_unfiltered.csv
    ├── rep2/<patient>_<timepoint>_rep2_soft.csv
    └── <patient>_<timepoint>_merged.csv
```

### 2. Depth Files Directory (iVar outputs)

```
ivar_runs/
└── <patient>/<timepoint>/<rep>/masked_variants_depth.tsv
```

### 3. Metadata File (required for plotting)

A `.csv` file with:

* `sample_name`: e.g. `N1_0_rep1_unfiltered`, or `N1-A` for acute samples
* `ct_value`: numeric cycle threshold
* Optional: `timepoint` (days since index case symptom onset)

---

## Usage

### 1. Run Diversity Analysis

```bash
python analyze_diversity.py \
    --mutation-analysis-path ./mutation_analysis \
    --ivar-runs-path ./ivar_runs \
    --ct-values-file ./metadata.csv \
    --output-path ./results \
    --project-name chronic
```

### 2. Output Files

| File                                                        | Description                                  |
| ----------------------------------------------------------- | -------------------------------------------- |
| `chronic_diversity_summary.csv`                             | Pivoted table of π and overdispersion values |
| `chronic_pi_diversity_plot.png`                             | Regression plot of π vs. Ct                  |
| `chronic_pi_diversity_plot_per_patient_grid_timecourse.png` | Timecourse π plots per patient               |
| `chronic_pi_diversity_plot_aggregated_timecourse.png`       | Aggregate π trends over time                 |
| `chronic_pi_diversity_plot_ct_over_time.png`                | Ct trajectories over time                    |

---

## Plot Details

* **Nucleotide Diversity (π)**: computed according to the formula introduced in Nei & Li, 1979 (DOI: 10.1073/pnas.76.10.5269).

Plots include:

* Regression line with 95% CI
* Per-patient longitudinal diversity
* Aggregate diversity across cohorts

