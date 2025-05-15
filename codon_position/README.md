# Codon Position Enrichment Analysis in SARS-CoV-2 Mutations

This script analyzes the distribution of mutations across codon positions (1st, 2nd, 3rd) in SARS-CoV-2 coding regions. It normalizes observed mutation counts by the reference genome's base composition and identifies statistically enriched codon positions using binomial tests.

- Loads filtered mutation data from iVar output (`masked_variants.tsv`)
- Maps mutations to codon positions based on GFF CDS annotations
- Normalizes mutation frequencies by base composition
- Computes statistical enrichment for each codon position
- Generates publication-quality plots and summary tables

## Usage

```bash
python analyze_codon_position_enrichment.py \
    --gff refs/wuhan.gff.gz \
    --fasta refs/wuhan.fasta \
    --input-dir chronic/ivar_runs/wuhan_ref_min_freq_0.0/ \
    --output-plot normalized_codon_position_bins.png
```

### Arguments

* `--gff`: Path to compressed GFF file with CDS annotations
* `--fasta`: Path to reference genome in FASTA format
* `--input-dir`: Directory containing `masked_variants.tsv` files (iVar output)
* `--output-plot`: Path to save the output PNG plot (raw data will be saved as a CSV alongside)

## Output

* A PNG plot showing normalized mutation proportions across codon positions for various frequency bins
* A CSV file with the raw data used in the plot (`<output-plot>_raw.csv`)
