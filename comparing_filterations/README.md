# Dropped Mutation Analysis

This script analyzes mutations filtered in the SARS-CoV-2 mutation pipeline and compares them with known **problematic** and **suspicious** genomic sites. It generates overlap visualizations and summarizes drop reasons.

> Requires input from the **original version** of the mutation filtering pipeline (before suspicious/problematic sites were added to the soft filter).

---

## Usage

```bash
python analyze_dropped_positions.py \
  --dropped-csv path/to/dropped_mutation_locations.csv \
  --suspicious-csv path/to/suspicious_mutations.csv \
  --output-dir path/to/output/
```

---

## Inputs

* `dropped_mutation_locations.csv`: Contains `mutation`, `drop_reason`, and `phase`
* `suspicious_mutations.csv`: CSV with a `mut_nuc` column (e.g. `"A123C"`)
* Problematic sites VCF: auto-downloaded

---

## Outputs

* `mutation_position_overlap_venn3.png`: Venn diagram
* `mutation_position_upset_4set.png`: UpSet plot
* Drop reason summaries: bar plots + CSVs for soft and replicate filters

---

Let me know if you'd like to generate these plots from multiple filtering runs or combine across cohorts.
