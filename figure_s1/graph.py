"""
Uses the result of prepare_mutation_input.py to create two graphs:
1. Proportions of mutations at different codon positions (y axis) by their frequency (x-axis) (used in figure S1)
2. Amount of mutations (y axis) by their frequency (x-axis)
"""
import click
import os
import pandas as pd
import plotly.graph_objects as go
import logging
from typing import List

from figure_s1.prepare_mutation_input import analyze_mutations, columns

logging.basicConfig(level=logging.DEBUG)

buckets = [1e-4, 0.5e-3, 1e-3, 0.5e-2, 1e-2, 0.5e-1]


def get_mutation_count_by_frequency(df: pd.DataFrame, output_directory: str) -> str:
    """
    Creates a graph of amount of mutations (y axis) by their frequency (x-axis)
    """
    count_fig = go.Figure()
    sorted_df = df.sort_values(by=columns.mutation_frequency)

    count_fig.add_trace(go.Scatter(
            x=sorted_df[columns.mutation_frequency],
            y=sorted_df['mutation_count_by_frequency'],
            mode='lines+markers',
            name=f'Mutation Count',
            marker=dict(symbol='circle')
    ))

    count_fig.update_layout(
        xaxis=dict(
            title='Mutation Frequency',
            type='log',
            dtick=1, 
            tickvals=buckets,
            tickformat='.0e',
        ),
        yaxis=dict(
            title='Mutation Count',
        ),
        title='Count of Mutation by Frequency of Mutation',
        template='plotly_white'
    )


    count_res_path = os.path.join(output_directory, 'mutation_count_by_frequency.html')
    logging.info(f'Saving plot under {count_res_path}')
    count_fig.write_html(count_res_path)
    return count_res_path

def get_mutation_proportion_per_codon_position_by_frequency(df: pd.DataFrame, output_directory: str) -> str:
    """
    Creates a graph of proportions of mutations at different codon positions (y axis) by their frequency (x-axis)
    """
    proportion_fig = go.Figure()

    # Plot each codon position
    for key, group in df.groupby(columns.mutated_position_in_codon):
        sorted_group = group.sort_values(by=columns.mutation_frequency)
        proportion_fig.add_trace(go.Scatter(
            x=sorted_group[columns.mutation_frequency],
            y=sorted_group['proportion'],
            mode='lines+markers',
            name=f'{key}',
            marker=dict(symbol='circle'),
            line=dict(width=4)
        ))
    
    proportion_fig.update_layout(
        xaxis=dict(
            title='Mutation Frequency',
            type='log',
            dtick=1, 
            tickvals=buckets,
            tickformat='.0e'
        ),
        yaxis=dict(
            title='Proportion',
            range=[0, 0.4]
        ),
        title='',
        legend_title='Position in Codon',
        template='plotly_white',
        font=dict(size=16)
    )

    proportion_res_path = os.path.join(output_directory, 'mutation_proportion_per_codon_position_by_frequency.html')
    logging.info(f'{proportion_res_path}: Saving plot under {proportion_res_path}')
    proportion_fig.write_html(proportion_res_path)
    return proportion_res_path

def get_plots(mutation_analysis_path: str, output_directory: str) -> List[str]:
    logging.info(f'{mutation_analysis_path}: Starting plot creation')

    df = pd.read_csv(mutation_analysis_path, usecols=[columns.mutated_position_in_codon, columns.mutation_frequency, columns.sequence_type])
    
    # filter coding sequence
    df = df[df[columns.sequence_type] == 'CDS']
    df.drop([columns.sequence_type], axis=1, inplace=True)

    logging.info(f"above 0.1: {len(df[df[columns.mutation_frequency] > 0.1])}")
    logging.info(f"between 0.1 and 0.01: {len(df[(df[columns.mutation_frequency] < 0.1) & (df[columns.mutation_frequency] > 0.01)])}")
    logging.info(f"below 0.01: {len(df[df[columns.mutation_frequency] < 0.01])}")
    
    df[columns.mutation_frequency] = pd.cut(df[columns.mutation_frequency], bins=buckets, labels=buckets[:-1])
    
    # Group by mutation frequency and codon position, then count occurrences
    mutation_count_by_frequency = df.groupby(columns.mutation_frequency).size().rename('mutation_count_by_frequency')
    codon_position_mutation_count_by_frequency = df.groupby([columns.mutation_frequency, columns.mutated_position_in_codon]).size().rename('mutation_count_by_frequency_and_codon_position')

    # Combine the counts into a single DataFrame
    merged_df = mutation_count_by_frequency.reset_index().merge(
        codon_position_mutation_count_by_frequency.reset_index(),
        on=columns.mutation_frequency
    )
    merged_df['proportion'] = merged_df['mutation_count_by_frequency_and_codon_position'] / merged_df['mutation_count_by_frequency']
    return [
        get_mutation_proportion_per_codon_position_by_frequency(merged_df, output_directory),
        get_mutation_count_by_frequency(merged_df, output_directory),
    ]

@click.command()
@click.option('-i', '--input-paths', type=click.Path(exists=True, readable=True), multiple=True, help='Directories containing input files to be analyzed. Files can either be freqs.tsv or masked.tsv')
@click.option('-csv', '--csv-input', type=click.Path(exists=True, readable=True, dir_okay=False), help='A CSV file describing mutations, generated by mutation_analysis.py. If provided, will be used for analysis instead of raw freqs files.')
@click.option('-o', '--output-directory', show_default=True, type=click.Path(exists=True, writable=True, dir_okay=True, file_okay=False), default=os.path.expanduser('~'), help='A directory in which to store the resulting CSV and plots')
def cli(input_paths: List[str], csv_input: str, output_directory: str) -> None:
    if not input_paths and not csv_input:
        raise click.UsageError("You must provide either input files or prepare_mutation_input output in csv form.")
    
    if input_paths and csv_input:
        raise click.UsageError("Please provide only one of input files or mutation analysis CSV file, not both.")

    if csv_input:
        logging.info(f'Loading mutation information from {csv_input}')
        mutation_analysis_path = csv_input

    if input_paths:
        logging.info(f'Preparing mutation input from {input_paths}')
        mutation_analysis_path = analyze_mutations(input_paths, output_directory)

    get_plots(mutation_analysis_path, output_directory)

if __name__ == '__main__':
    cli()
