# 🧬 Suspicious SARS-CoV-2 Mutations Pipeline

This repository provides a full Python pipeline for identifying and analyzing **suspicious recurrent low-frequency mutations** in SARS-CoV-2 deep sequencing data. It includes variant filtering, recurrence detection, annotation, primer distance analysis, gene enrichment, and statistical testing.

---

## 📁 Overview

### Main Features:
- 📦 Combines masked variant files from iVar output folders
- 🔎 Filters mutations based on frequency and coverage
- 📊 Identifies mutations enriched within VOCs (Variants of Concern)
- 🧬 Annotates mutations with ORF and amino acid changes using reference genome
- 🔬 Calculates proximity to V3/V4 primers
- 🧪 Performs gene-level enrichment (hypergeometric test)
- 🧠 Visualizes mutation fitness and primer distance distributions
- 📈 Applies Mann–Whitney U test comparing suspicious vs. background distances

---

## 📂 Input Files

| Type                | Description                                                       |
|---------------------|-------------------------------------------------------------------|
| `ivar_trim_X/`      | Base folders containing masked variant TSVs per sample            |
| `sra_metadata.csv`  | Metadata file with `srr`, `voc`, and `sample` columns             |
| `V3_nCoV_2019.tsv`  | Primer positions (start, end, strand) for V3 primer set           |
| `V4_nCoV_2019.tsv`  | Same for V4 primer set                                            |
| `wuhan_ref.fasta`   | Wuhan reference genome for annotation                             |
| `fitness.tsv`       | Mutation-level fitness scores and synonymity annotations          |


---

## 🚀 Usage

Run the script directly:

```bash
python suspicious_mutation_analysis.py
```

All output will be printed and saved in the current directory.

You may update `base_paths`, `metadata_path`, and other paths inside the `if __name__ == "__main__"` block.

---

## 📤 Outputs

- `suspected_muts_with_soft_filtering_trim_*.csv`: Annotated suspicious mutations
- Enrichment table (printed to console)
- Interactive Plotly figure (opened in browser)
- Mann–Whitney U test results (printed for V3 and V4)

---

## 📊 Methods Summary

### 🧼 Filtering Criteria
- Frequency: 1% ≤ ALT_FREQ ≤ 20%
- Depth: ALT_DP ≥ 50, TOTAL_DP > 100

### 🔁 Recurrence Logic
- A mutation is suspicious if found in ≥ 20% of samples within a VOC

### 🧠 Annotation
- Annotated using codon-based logic with Wuhan reference genome
- Supports SNPs, insertions, and deletions

### 🔬 Enrichment
- Hypergeometric test comparing number of mutations per gene to genome-wide expectation

### 📈 Visualization
- KDE plots of fitness and primer distance (mutations vs. genome-wide background)

---

## 📌 Example Output

```bash
--- Processing: trim_40 ---
X suspicious mutations saved for trim_40

Gene enrichment:
╒══════════╤════╤════╤══════╤══════╤══════════╤══════════════════╤════════════════════╕
│ Gene     │ k  │ s  │ M    │ N    │ Expected │ Fold Enrichment  │ Hypergeometric p   │
├──────────┼────┼────┼──────┼──────┼──────────┼──────────────────┼────────────────────┤
│ S        │    │    │      │      │          │                  │                    │
│ N        │    │    │      │      │          │                  │                    │
│ ORF1ab   │    │    │      │      │          │                  │                    │
╘══════════╧════╧════╧══════╧══════╧══════════╧══════════════════╧════════════════════╛

Density plot saved to export_path/my_density_plot.png

Mann-Whitney test results:
V3: {'U statistic': 12345.0, 'p-value': 0.00018, ...}
V4: {'U statistic': 12345.0, 'p-value': 0.00429, ...}
```

---

## 📄 License

MIT License  
Feel free to reuse or modify.
