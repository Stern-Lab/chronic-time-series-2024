"""
Arranges information about mutations in SARS-CoV-2 sequences - adding information about codon position and protein
Input:
    #    Directories containing input frequency files. Files can either be from PyAccuNGS pipeline (freqs.tsv) or from iVar pipeline (masked.tsv).
Ouput:
    #    A CSV file containing all mutations found with their frequency, adding data about their position in the codon and protein.
"""
import click
import os
import numpy as np
import pandas as pd
import logging

from Bio.Seq import Seq
from typing import Dict, List, Tuple, Optional
import dask.dataframe as dd
from .utils import FREQS_FILE_SUFFIX, IVAR_FILE_SUFFIX, MutationAnalysisColumns, convert_codon_dict_to_string, get_reference_codon, find_frequency_files

logging.basicConfig(level=logging.INFO)
columns = MutationAnalysisColumns()

@click.command()
@click.option('-i', '--input', 'input_directories', type=click.Path(exists=True, readable=True), multiple=True, help='Directories containing input files to be analyzed. Files can either be freqs.tsv or masked.tsv')
@click.option('-o', '--output', 'output_directory', type=click.Path(exists=True, writable=True, file_okay=False, dir_okay=True), default=os.path.join(os.path.expanduser('~')), help="Path to save the aggregated mutation data.")
def main(input_directories: List[str], output_directory: str):
    analyze_mutations(input_directories, output_directory)

def analyze_mutations(input_directories: List[str], output_directory: str) -> str:
    files = find_frequency_files(input_directories)
    output_file = create_mutation_analysis_multi_sample(files, output_directory)
    logging.info(f"Aggregated mutation data saved to {output_file}")
    return output_file

def create_mutation_analysis_multi_sample(files: List[str], output_directory: str) -> str:
    """Performs mutation analysis across multiple samples and saves to a CSV file."""
    logging.info(f"Performing analysis on {len(files)} samples")

    # Create a Dask DataFrame using from_map
    combined_df = dd.from_map(
        lambda file: _create_mutation_analysis_single_sample(file, output_directory),
        files,
        meta=columns.meta
    ).compute()

    output_file_path = os.path.join(output_directory, 'mutation_analysis.csv')
    logging.info(f'Saving mutation analysis in {output_file_path}')
    combined_df.to_csv(output_file_path, index=False)
    
    logging.info(f'Finished mutation analysis')
    return output_file_path

def _process_mutated_codon(row: pd.Series, reference_codon: Dict[int, str], codon_positions: Tuple[int, int, int]) -> Optional[pd.Series]:
    """Processes a single row of mutation data and returns a pandas Series."""
    protein_data = columns.get_protein_data(codon_positions[0])
    mutated_codon_dict = reference_codon.copy()
    mutated_position = row['ref_pos']
    mutated_position_in_codon = mutated_position % 3 or 3
    mutated_base = row['read_base']
    frequency = row['frequency']

    mutated_codon_dict[mutated_position] = mutated_base
    mutated_codon_string = convert_codon_dict_to_string(mutated_codon_dict)
    reference_codon_string = convert_codon_dict_to_string(reference_codon)
    reference_amino_acid = Seq(reference_codon_string).translate()
    mutated_amino_acid = Seq(mutated_codon_string).translate()

    return pd.Series({
        columns.mutated_position: mutated_position,
        columns.mutation_string: reference_amino_acid + str(mutated_position) + mutated_amino_acid,
        columns.reference_codon: reference_codon_string,
        columns.mutated_codon: mutated_codon_string,
        columns.mutated_position_in_codon: mutated_position_in_codon,
        columns.mutation_frequency: frequency,
        columns.reference_amino_acid: reference_amino_acid,
        columns.mutated_amino_acid: mutated_amino_acid,
        **protein_data
    })

def _process_codon(codon_positions: Tuple[int, int, int], freqs_data: pd.DataFrame) -> dd.DataFrame:
    """Processes mutations for a given codon."""
    codon_data = freqs_data[freqs_data['ref_pos'].isin(codon_positions)]
    reference_codon = get_reference_codon(codon_positions, codon_data)
    mutation_data = codon_data[codon_data['ref_base'] != codon_data['read_base']]
    columns_order = columns.all_columns
    columns_order.remove(columns.sample)
    columns_order.remove(columns.file_type)
    return mutation_data.apply(lambda row: _process_mutated_codon(row, reference_codon, codon_positions), axis=1).reindex(columns=columns_order)

def _process_freqs(file_path: str) -> dd.DataFrame:
    """
    Processes PyAccuNGS files.
    In these files we don't have the codon information, only mutation frequency per position.
    Therefore, we create a set of all codons in the sequences positions, and arrange the mutation information by codon.
    """
    sample_number = os.path.basename(os.path.dirname(file_path))
    df = pd.read_csv(file_path, sep='\t', usecols=['ref_pos', 'ref_base', 'read_base', 'frequency'])
    df = df[df['ref_pos'].apply(lambda x: x.is_integer())]  # Remove insertions
    df['ref_pos'] = df['ref_pos'].astype(int)
    df = df[df['read_base'] != '-']  # Remove deletions
    df = df[df['frequency'] > 0]
    logging.debug(f"{file_path}: After filter, got {len(df)} rows")
    ref_pos_set = set(df['ref_pos'])
    codons = [
        (start_position, start_position + 1, start_position + 2)
        for start_position in range(1, int(max(df['ref_pos'])) - 1, 3)
        if {start_position, start_position + 1, start_position + 2}.issubset(ref_pos_set)
    ]

    columns_meta = columns.meta
    columns_meta.pop(columns.sample)
    columns_meta.pop(columns.file_type)

    # Process each codon in parallel
    result = dd.from_map(lambda codon: _process_codon(codon, df), codons, meta=columns_meta).compute()
    logging.debug(f"{file_path}: Got {len(result)} mutations")

    result[columns.sample] = sample_number
    result[columns.file_type] = FREQS_FILE_SUFFIX
    result = result.reindex(columns=columns.all_columns)
    return result

def _process_ivar(file_path: str) -> pd.DataFrame:
    """
    Processes iVar frequency files.
    In these files we already have the matching codon for each position, so we just add the protein information.
    """
    sample_number = os.path.basename(file_path)[:-(len(IVAR_FILE_SUFFIX)+1)]
    df = pd.read_csv(file_path, sep='\t', usecols=['POS', 'REF', 'ALT', 'ALT_FREQ', 'REF_CODON', 'ALT_CODON', 'REF_AA', 'ALT_AA'])
    df = df[df['ALT'].apply(lambda x: len(x) == 1 and '+' not in x and '-' not in x)]  # Remove insertions and deletions
    df[columns.mutation_string] = df['REF_AA'] + df['POS'].astype(str) + df['ALT_AA']

    df.drop(['REF', 'ALT'], axis=1, inplace=True)
    df.rename(
        columns={
            'POS': columns.mutated_position,
            'ALT_FREQ': columns.mutation_frequency,
            'REF_CODON': columns.reference_codon,
            'ALT_CODON': columns.mutated_codon,
            'REF_AA': columns.reference_amino_acid,
            'ALT_AA': columns.mutated_amino_acid
        },
        inplace=True
    )

    df[columns.mutated_position] = df[columns.mutated_position].astype(int)
    df[columns.mutated_position_in_codon] = np.where(df[columns.mutated_position] % 3 == 0, 3, df[columns.mutated_position] % 3)
    df[columns.mutation_frequency] = df[columns.mutation_frequency].astype(float)

    # Add protein information
    protein_info = df[columns.mutated_position].apply(columns.get_protein_data)
    protein_info_df = pd.DataFrame(protein_info.tolist(), index=protein_info.index)
    df = pd.concat([df, protein_info_df], axis=1)
    
    df[columns.sample] = sample_number
    df[columns.file_type] = IVAR_FILE_SUFFIX
    
    # Ensure the column types match the expected metadata types
    for col, dtype in columns.meta.items():
        df[col] = df[col].astype(dtype)
    
    df = df.reindex(columns=columns.all_columns)
    return df

def _create_mutation_analysis_single_sample(file_path: str, output_directory: str):
    logging.info(f"Processing {file_path}")
    if file_path.endswith(FREQS_FILE_SUFFIX):
        result = _process_freqs(file_path)
    elif file_path.endswith(IVAR_FILE_SUFFIX):
        result = _process_ivar(file_path)
    else:
        raise ValueError(f"Unrecognized file type: {file_path}")
    
    output_file = os.path.join(output_directory, f"{result[columns.sample].iloc[0]}_mutation_analysis.csv")
    result.to_csv(output_file, index=False)
    logging.info(f"Finished processing {file_path}")
    return result


if __name__ == "__main__":
    main()