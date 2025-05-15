# Nucleotide Diversity (π) Comparison Plot

This script visualizes SARS-CoV-2 nucleotide diversity (`π`) in acute vs. chronic samples after replicate concordance filtering. It performs a statistical test and outputs a log-scaled boxplot.

- Reads summary CSVs from acute and chronic sample analyses
- Computes a Mann–Whitney U test on `pi_merged` values
- Outputs a publication-quality boxplot (log-scaled y-axis, annotated p-value)

## Usage

```bash
python plot_pi_diversity_boxplot.py \
    --acute /path/to/acute_diversity_summary.csv \
    --chronic /path/to/chronic_diversity_summary.csv \
    --output nucleotide_diversity_boxplot.png
```

### Arguments

* `--acute`: Path to CSV file containing acute diversity data (must include a `pi_merged` column)
* `--chronic`: Path to CSV file containing chronic diversity data (same format)
* `--output`: Filename for the output PNG (default: `nucleotide_diversity_boxplot.png`)

## Notes

* `pi_merged` should reflect replicate-concordance-filtered values from your diversity pipeline.