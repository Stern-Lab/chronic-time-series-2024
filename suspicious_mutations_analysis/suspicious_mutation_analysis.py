import os
import re
import pandas as pd
import numpy as np
from typing import List, Tuple, Optional, Dict
from Bio import SeqIO, Seq
from tabulate import tabulate
from scipy.stats import hypergeom, mannwhitneyu
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde

# === CONSTANTS ===
ORF_COORDS = {
    'ORF1a': (266, 13483),
    'ORF1b': (13468, 21555),
    'S': (21563, 25384),
    'ORF3a': (25393, 26220),
    'E': (26245, 26472),
    'M': (26523, 27191),
    'ORF6': (27202, 27387),
    'ORF7a': (27394, 27759),
    'ORF7b': (27756, 27887),
    'ORF8': (27894, 28259),
    'N': (28274, 29533),
    'ORF10': (29558, 29674)
}
## Part 2: Data Loading & Filtering
def combine_masked_variants(base_path: str) -> pd.DataFrame:
    """Combines variant files from subdirectories under a given base path."""
    all_dfs = []
    for subfolder in os.listdir(base_path):
        subfolder_path = os.path.join(base_path, subfolder)
        masked_variants_path = os.path.join(subfolder_path, 'masked_variants')
        if os.path.isdir(masked_variants_path):
            tsv_files = [f for f in os.listdir(masked_variants_path) if f.endswith('.tsv')]
            if tsv_files:
                try:
                    df = pd.read_csv(os.path.join(masked_variants_path, tsv_files[0]), sep='\t')
                    df.insert(0, 'srr', subfolder)
                    all_dfs.append(df)
                except pd.errors.ParserError as e:
                    print(f"Error reading {tsv_files[0]}: {e}")
    return pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()

def get_masked_trim_by_required_rules(
    df: pd.DataFrame,
    min_freq: float = 0.01,
    max_freq: float = 0.2,
    min_total_depth: int = 100,
    min_coverage_of_mut: int = 10
) -> pd.DataFrame:
    """Filters variants based on frequency and coverage thresholds."""
    filtered = df[
        (df['ALT_FREQ'] >= min_freq) &
        (df['ALT_FREQ'] <= max_freq) &
        (df['TOTAL_DP'] > min_total_depth) &
        (df['ALT_DP'] > min_coverage_of_mut)
    ].copy()

    filtered = filtered[['srr', 'POS', 'REF', 'ALT', 'REF_DP', 'REF_RV', 'REF_QUAL',
                         'ALT_DP', 'ALT_RV', 'ALT_QUAL', 'ALT_FREQ', 'TOTAL_DP', 'PVAL', 'PASS']]
    filtered['mut_nuc'] = filtered['REF'] + filtered['POS'].astype(str) + filtered['ALT']
    return filtered

def get_common_mutations(df: pd.DataFrame, threshold: float, sample_size: int) -> List[str]:
    """Returns mutations that appear in at least a threshold fraction of samples."""
    mutation_counts = df.groupby('mut_nuc')['srr'].nunique()
    mutation_freq = mutation_counts / sample_size
    return mutation_freq[mutation_freq >= threshold].index.tolist()

# Mutation Summaries & Annotation
def calculate_averages_filtered(df: pd.DataFrame, suspicious_mut_list: List[str]) -> Tuple[float, float, float]:
    """Calculates mean ALT_DP, TOTAL_DP, and ALT_FREQ for suspicious mutations."""
    df_filtered = df[df['mut_nuc'].isin(suspicious_mut_list)]
    if df_filtered.empty:
        return 0.0, 0.0, 0.0
    return (
        df_filtered['ALT_DP'].mean(),
        df_filtered['TOTAL_DP'].mean(),
        df_filtered['ALT_FREQ'].mean()
    )

def extract_pos(mutation: str) -> Optional[int]:
    match = re.search(r'\d+', str(mutation))
    return int(match.group()) if match else None

def annotate_mutation(mutation: str, ref_seq: str) -> Tuple[Optional[str], str, str]:
    match_sub = re.match(r'^([ACGT])(\d+)([ACGT])$', mutation)
    match_ins = re.match(r'^([ACGT])(\d+)\+([ACGT]+)$', mutation)
    match_del = re.match(r'^([ACGT])(\d+)-([ACGT]*)$', mutation)

    if match_sub:
        ref, pos, alt = match_sub.groups()
        pos = int(pos)
        if ref_seq[pos - 1].upper() != ref:
            return None, mutation, 'mismatch'
        for orf, (start, end) in ORF_COORDS.items():
            if start <= pos <= end:
                rel_pos = pos - start
                codon_start = start + (rel_pos // 3) * 3
                codon_seq = ref_seq[codon_start - 1: codon_start + 2]
                if len(codon_seq) != 3:
                    return orf, mutation, 'syn'
                codon_list = list(codon_seq)
                codon_pos = pos - codon_start
                codon_list[codon_pos] = alt
                orig_codon = Seq.Seq(codon_seq)
                new_codon = Seq.Seq(''.join(codon_list))
                aa_orig = orig_codon.translate()
                aa_new = new_codon.translate()
                aa_pos = (rel_pos // 3) + 1
                return (orf, f"{aa_orig}{aa_pos}{aa_new}", 'syn' if aa_orig == aa_new else 'non-syn')
        return None, mutation, 'syn'

    elif match_ins or match_del:
        pos = int(match_ins.groups()[1] if match_ins else match_del.groups()[1])
        for orf, (start, end) in ORF_COORDS.items():
            if start <= pos <= end:
                return orf, mutation, 'indel'
        return None, mutation, 'indel'

    return None, mutation, 'unknown'

def annotate_dataframe(df: pd.DataFrame, ref_seq: str) -> pd.DataFrame:
    results = df['mut_nuc'].apply(lambda x: annotate_mutation(x, ref_seq))
    df[['orf', 'mut_AA', 'mut_kind']] = pd.DataFrame(results.tolist(), index=df.index)
    return df

def load_wuhan_ref(path_to_fasta: str) -> str:
    record = SeqIO.read(path_to_fasta, "fasta")
    return str(record.seq).upper()
def get_variant_specific_common_mutations(
    df: pd.DataFrame,
    threshold: float,
    mdf_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Identifies mutations present in ≥ threshold fraction of samples per VOC,
    and summarizes ALT_DP, TOTAL_DP, ALT_FREQ per mutation.

    Args:
        df: DataFrame with columns ['mut_nuc', 'srr', 'voc', 'ALT_DP', 'TOTAL_DP', 'ALT_FREQ'].
        threshold: Proportion threshold (0 < threshold ≤ 1).
        mdf_df: Metadata DataFrame with ['srr', 'voc'] to determine total samples per VOC.

    Returns:
        pd.DataFrame: Summary of suspicious mutations including per-VOC counts and averages.
    """
    vocs = df['voc'].dropna().unique()
    
    # Total number of unique SRRs per VOC
    num_patients_per_voc = (
        mdf_df.groupby('voc')['srr']
              .nunique()
              .to_dict()
    )

    suspected_mutations = set()

    for voc in vocs:
        df_voc = df[df['voc'] == voc]
        num_patients = num_patients_per_voc.get(voc, 0)
        if num_patients == 0:
            continue
        mut_counts = df_voc.groupby('mut_nuc')['srr'].nunique()
        mut_freq = mut_counts / num_patients
        passing_muts = mut_freq[mut_freq >= threshold].index.tolist()
        suspected_mutations.update(passing_muts)

    # Filter and summarize
    filtered_df = df[df['mut_nuc'].isin(suspected_mutations)].copy()
    summary = pd.DataFrame({'mut_nuc': list(suspected_mutations)})

    for voc in vocs:
        voc_counts = (
            filtered_df[filtered_df['voc'] == voc]
            .groupby('mut_nuc')
            .size()
            .reindex(summary['mut_nuc'], fill_value=0)
            .reset_index(drop=True)
        )
        summary[f"{voc}_mut_count"] = voc_counts

    # Add average depth and frequency
    avg_metrics = (
        filtered_df.groupby('mut_nuc')[['ALT_DP', 'TOTAL_DP', 'ALT_FREQ']]
        .mean()
        .reindex(summary['mut_nuc'])
        .reset_index(drop=True)
    )

    summary = pd.concat([summary, avg_metrics], axis=1)
    return summary

# Distance to Primers & Fitness Merge
def calculate_min_distances_with_fitness(
    mutation_df: pd.DataFrame,
    v3_df: pd.DataFrame,
    v4_df: pd.DataFrame,
    metadata_df: pd.DataFrame
) -> pd.DataFrame:
    def extract_numeric_position(mut_str: str) -> Optional[int]:
        digits = re.findall(r'\d+', str(mut_str))
        return int(digits[0]) if digits else None

    def get_min_distance_and_strand(mut_pos: int, primers_df: pd.DataFrame) -> Tuple[int, str]:
        distances = primers_df.apply(
            lambda row: {
                'distance': min(abs(mut_pos - row['start']), abs(mut_pos - row['end'])),
                'strand': row['strand']
            }, axis=1
        )
        min_result = min(distances, key=lambda x: x['distance'])
        return min_result['distance'], min_result['strand']

    mutation_df['mut_pos'] = mutation_df['mut_nuc'].apply(extract_numeric_position)

    v3_results = mutation_df['mut_pos'].apply(lambda pos: get_min_distance_and_strand(pos, v3_df) if pd.notnull(pos) else (None, None))
    v4_results = mutation_df['mut_pos'].apply(lambda pos: get_min_distance_and_strand(pos, v4_df) if pd.notnull(pos) else (None, None))

    mutation_df[['V3_distance', 'V3_strand']] = pd.DataFrame(v3_results.tolist(), index=mutation_df.index)
    mutation_df[['V4_distance', 'V4_strand']] = pd.DataFrame(v4_results.tolist(), index=mutation_df.index)

    metadata_df = metadata_df.rename(columns={'mutation': 'mut_nuc'})
    mutation_df = mutation_df.merge(
        metadata_df[['mut_nuc', 'delta_fitness', 'synonymity']],
        how='left', on='mut_nuc'
    )
    return mutation_df

# Enrichment Analysis
def enrichment_row(gene: str, k: int, s: int, M: int, N: int) -> List:
    expected = (M / N) * s
    fold = k / expected if expected != 0 else float('inf')
    pval = hypergeom.sf(k - 1, N, M, s)
    return [gene, k, s, M, N, round(expected, 2), round(fold, 2), f"{pval:.2e}"]

def perform_gene_enrichment(filtered_df: pd.DataFrame, genome_size: int = 29903) -> pd.DataFrame:
    gene_lengths = {
        "S": 25384 - 21563 + 1,
        "N": 29533 - 28274 + 1,
        "ORF1a": 13468 - 266 + 1,
        "ORF1b": 21555 - 13468 + 1
    }
    all_coding = sum(gene_lengths.values())
    s = len(filtered_df)
    k_s = filtered_df[filtered_df['orf'] == 'S'].shape[0]
    k_n = filtered_df[filtered_df['orf'] == 'N'].shape[0]
    k_orf1ab = filtered_df[filtered_df['orf'].isin(['ORF1a', 'ORF1b'])].shape[0]

    table = [
        enrichment_row("S", k_s, s, gene_lengths['S'], genome_size),
        enrichment_row("N", k_n, s, gene_lengths['N'], genome_size),
        enrichment_row("ORF1ab", k_orf1ab, s, gene_lengths['ORF1a'] + gene_lengths['ORF1b'], genome_size)
    ]

    headers = ["Gene", "k", "s", "M", "N", "Expected", "Fold Enrichment", "Hypergeometric p"]
    return pd.DataFrame(table, columns=headers)

# Density Plotting

def compute_genome_distances_fast(ref_seq: str, primers_df: pd.DataFrame) -> pd.DataFrame:
    genome_positions = np.arange(1, len(ref_seq) + 1)
    dists_to_starts = np.abs(genome_positions[:, None] - primers_df['start'].values[None, :])
    dists_to_ends = np.abs(genome_positions[:, None] - primers_df['end'].values[None, :])
    minimal_dists = np.minimum(dists_to_starts, dists_to_ends).min(axis=1)
    on_primer = np.any(
        (primers_df['start'].values[None, :] <= genome_positions[:, None]) &
        (primers_df['end'].values[None, :] >= genome_positions[:, None]),
        axis=1
    )
    return pd.DataFrame({
        'pos': genome_positions,
        'minimal_dist': minimal_dists,
        'on_primer': on_primer
    })

def mann_whitney_test(bases_df: pd.DataFrame, mutations_df: pd.DataFrame) -> Dict[str, float]:
    """Performs a one-sided Mann-Whitney U test (mutations > background)."""
    base_dists = bases_df[bases_df['on_primer'] == False]['minimal_dist']
    mut_dists = mutations_df[mutations_df['on_primer'] == False]['minimal_dist']
    u, p = mannwhitneyu(base_dists, mut_dists, alternative='greater')
    return {
        'U statistic': u,
        'p-value': p,
        'control_mean': base_dists.mean(),
        'control_std': base_dists.std(),
        'mutations_mean': mut_dists.mean(),
        'mutations_std': mut_dists.std()
    }
def plot_full_density_figure(
    filtered_df: pd.DataFrame,
    ref_seq: str,
    v3_df: pd.DataFrame,
    v4_df: pd.DataFrame,
    width: int = 1200,
    height: int = 700
) -> go.Figure:
    """
    Creates a 2x3 Plotly figure comparing:
    - Fitness density by mutation kind
    - Distance to V3/V4 primers (mutations vs. genome-wide)

    Args:
        filtered_df: Annotated DataFrame with suspicious mutations
        ref_seq: Full reference sequence as string
        v3_df: Primer table for V3 (with 'start', 'end', 'strand')
        v4_df: Primer table for V4
        width: Plot width
        height: Plot height

    Returns:
        Plotly Figure object
    """

    def compute_genome_distances(ref_seq: str, primers_df: pd.DataFrame) -> pd.DataFrame:
        genome_positions = np.arange(1, len(ref_seq) + 1)
        dists_to_starts = np.abs(genome_positions[:, None] - primers_df['start'].values[None, :])
        dists_to_ends = np.abs(genome_positions[:, None] - primers_df['end'].values[None, :])
        minimal_dists = np.minimum(dists_to_starts, dists_to_ends).min(axis=1)
        on_primer = np.any(
            (primers_df['start'].values[None, :] <= genome_positions[:, None]) &
            (primers_df['end'].values[None, :] >= genome_positions[:, None]),
            axis=1
        )
        return pd.DataFrame({'pos': genome_positions, 'minimal_dist': minimal_dists, 'on_primer': on_primer})

    def density_trace(data: pd.Series, name: str, line_color: str, fill_color: str, opacity: float = 1.0, showlegend: bool = True) -> go.Scatter:
        kde = gaussian_kde(data)
        x_vals = np.linspace(data.min() - 10, data.max() + 10, 500)
        y_vals = kde(x_vals)
        return go.Scatter(
            x=x_vals, y=y_vals, mode='lines', name=name,
            line=dict(color=line_color), opacity=opacity,
            fill='tozeroy', fillcolor=fill_color if fill_color else line_color,
            showlegend=showlegend
        )

    # === Compute genome-wide primer distances ===
    genome_v3_df = compute_genome_distances(ref_seq, v3_df)
    genome_v4_df = compute_genome_distances(ref_seq, v4_df)

    # === Subset mutation data ===
    fit_data = filtered_df.dropna(subset=['delta_fitness'])

    fig = make_subplots(
        rows=2, cols=3,
        column_widths=[0.4, 0.3, 0.3],
        row_heights=[0.55, 0.45],
        vertical_spacing=0.08,
        specs=[[{"rowspan": 2}, {}, {}], [None, {}, {}]],
        subplot_titles=[
            "Fitness (by Mutation Kind)",
            "Suspicious Mutations (V3)", "Genome-wide (V3)",
            "", "Suspicious Mutations (V4)", "Genome-wide (V4)"
        ]
    )

    # === Plot fitness by kind ===
    color_map = {
        'non-syn': ('blue', 'lightblue'),
        'syn': ('green', 'lightgreen'),
        'indel': ('deeppink', 'pink')
    }

    for mut_type, (line_c, fill_c) in color_map.items():
        subset = fit_data[fit_data['mut_kind'] == mut_type]['delta_fitness']
        if not subset.empty:
            fig.add_trace(density_trace(subset, mut_type, line_c, fill_c, opacity=0.4), row=1, col=1)

    # === Suspicious mutation distances (V3/V4) ===
    for i, col in enumerate(['V3_distance', 'V4_distance']):
        data = filtered_df[col].dropna()
        fig.add_trace(density_trace(data, 'Suspicious Mutations', 'purple', 'violet'), row=1 + (i > 0), col=2)

    # === Genome-wide distances (V3/V4) ===
    for i, genome_df in enumerate([genome_v3_df, genome_v4_df]):
        data = genome_df[genome_df['on_primer'] == False]['minimal_dist']
        fig.add_trace(density_trace(data, 'Genome-wide', 'steelblue', 'lightsteelblue', opacity=0.5), row=1 + (i > 0), col=3)

    # === Axes and layout ===
    fig.update_xaxes(title_text='Δ Fitness', row=1, col=1)
    fig.update_xaxes(title_text='Distance (nt)', row=1, col=2)
    fig.update_xaxes(title_text='Distance (nt)', row=2, col=2)
    fig.update_xaxes(title_text='Distance (nt)', row=1, col=3)
    fig.update_xaxes(title_text='Distance (nt)', row=2, col=3)

    fig.update_layout(
        height=height,
        width=width,
        template="plotly_white",
        showlegend=True,
        margin=dict(l=30, r=30, t=40, b=30),
        font=dict(size=12)
    )

    return fig

if __name__ == "__main__":
    # === Define paths ===
    base_paths = {
        "trim_20": "path/to/ivar_trim_20",
        "trim_30": "path/to/ivar_trim_30",
        "trim_40": "path/to/ivar_trim_40",
    }
    metadata_path = "path/to/sra_metadata.csv"
    ref_fasta = "path/to/wuhan_ref.fasta"
    primers_v3_path = "path/to/V3_nCoV_2019.tsv"
    primers_v4_path = "path/to/V4_nCoV_2019.tsv"
    fitness_path = "path/to/fitness.tsv"
    export_path = "path/to/export/figures"

    # === Load inputs ===
    md_df = pd.read_csv(metadata_path)
    ref_seq = load_wuhan_ref(ref_fasta)
    v3_df = pd.read_csv(primers_v3_path, sep="\t")
    v4_df = pd.read_csv(primers_v4_path, sep="\t")
    fitness_df = pd.read_csv(fitness_path, sep="\t")

    # === Process each trim level ===
    for trim_level, path in base_paths.items():
        print(f"\n--- Processing: {trim_level} ---")

        raw_df = combine_masked_variants(path)
        strict_filtered = get_masked_trim_by_required_rules(raw_df, 0.01, 0.2, 100, 50)
        strict_filtered = strict_filtered.merge(md_df[['srr', 'voc', 'primers_set']], on='srr', how='left')

        suspicious_df = get_variant_specific_common_mutations(strict_filtered, 0.2, md_df)
        annotated_df = annotate_dataframe(suspicious_df, ref_seq)
        full_info_df = calculate_min_distances_with_fitness(annotated_df, v3_df, v4_df, fitness_df)

        # Save
        full_info_df.to_csv(f"suspected_muts_with_soft_filtering_{trim_level}.csv", index=False)
        print(f"{len(full_info_df)} suspicious mutations saved for {trim_level}")

        # === Gene Enrichment ===
        enrichment_df = perform_gene_enrichment(full_info_df)
        print("\nGene enrichment:")
        print(tabulate(enrichment_df, headers='keys', tablefmt='fancy_grid'))

        # === Plotting ===
        fig = plot_full_density_figure(full_info_df, ref_seq, v3_df, v4_df, width=1200, height=700)
        fig.show()
        fig.write_image(f"{export_path}/my_density_plot.png")
        print(f"Density plot saved to {export_path}/my_density_plot.png")

        # === Mann-Whitney U ===
        genome_v3_df = compute_genome_distances_fast(ref_seq, v3_df)
        genome_v4_df = compute_genome_distances_fast(ref_seq, v4_df)

        mut_v3 = full_info_df[['V3_distance']].copy()
        mut_v3.columns = ['minimal_dist']
        mut_v3['on_primer'] = False

        mut_v4 = full_info_df[['V4_distance']].copy()
        mut_v4.columns = ['minimal_dist']
        mut_v4['on_primer'] = False

        print("\nMann-Whitney test results:")
        v3_stats = mann_whitney_test(genome_v3_df, mut_v3)
        v4_stats = mann_whitney_test(genome_v4_df, mut_v4)

        print("V3:", v3_stats)
        print("V4:", v4_stats)

