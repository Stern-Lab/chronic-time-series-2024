from dataclasses import dataclass
import io
import logging
import os
import threading
from typing import Dict, List, Tuple
from Bio import SeqIO
import requests

FREQS_FILE_SUFFIX = "freqs.tsv"
IVAR_FILE_SUFFIX = "masked.tsv"
# SARS-CoV-2 reference genome (NC_045512.2) in GenBank format
SARS_COV_2_NCBI = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&id=NC_045512.2&report=genbank&from=begin&to=end"

class GenBankDataSingleton:
    """A utility class for fetching and reading genbank data efficiently"""
    _instance = None
    _genbank_data = None
    _lock = threading.Lock()

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                # Double check if another thread has not already created an instance
                if cls._instance is None:
                    logging.info("Fetching GenBank data from NCBI...")
                    response = requests.get(SARS_COV_2_NCBI)
                    data = io.StringIO(response.text)
                    cls._instance = super(GenBankDataSingleton, cls).__new__(cls)
                    cls._instance._genbank_data = SeqIO.read(data, "genbank")
        return cls._instance
    
    def read(self) -> SeqIO.SeqRecord:
        with self._lock:
            return self._genbank_data
        
@dataclass
class MutationAnalysisColumns:
    """A utility class describing the columns of the resulting dataframe/csv of the mutation analysis"""
    mutated_position: str = "mutated_position"
    mutated_position_in_codon: str = "mutated_position_in_codon"
    mutation_string: str = "mutation_string"
    reference_codon: str = "reference_codon"
    mutated_codon: str = "mutated_codon"
    reference_amino_acid: str = "reference_amino_acid"
    mutated_amino_acid: str = "mutated_amino_acid"
    sequence_type: str = "sequence"
    protein_id: str = "protein_id"
    product: str = "product"
    gene: str = "gene"
    sample: str = "sample"
    mutation_frequency: str = "mutation_frequency"
    file_type: str = "file_type"

    @property
    def all_columns(self) -> List[str]:
        return list(self.meta.keys())

    @property
    def meta(self) -> Dict[str, str]:
        return {
            self.mutated_position: 'int64',
            self.mutated_position_in_codon: 'int64',
            self.mutation_string: 'string',
            self.reference_codon: 'string',
            self.mutated_codon: 'string',
            self.reference_amino_acid: 'string',
            self.mutated_amino_acid: 'string',
            self.sequence_type: 'string',
            self.protein_id: "string",
            self.product: "string",
            self.gene: "string",
            self.sample: 'string',
            self.mutation_frequency: 'float64',
            self.file_type: 'string'
        }
    
    def get_protein_data(self, position: int) -> Dict[str, str]:
        """Finds the protein information for a given nucleotide position."""
        unknown_string = 'Unknown'
        result = {
            self.sequence_type: unknown_string, 
            self.protein_id: unknown_string, 
            self.product: unknown_string, 
            self.gene: unknown_string
        }
        for feature in GenBankDataSingleton().read().features:
            if feature.location.start <= position <= feature.location.end:
                result[self.sequence_type] = feature.type
                if feature.type == "CDS":
                    result[self.protein_id] = feature.qualifiers.get('protein_id', [unknown_string])[0]
                    result[self.gene] = feature.qualifiers.get('gene', [unknown_string])[0]
                    result[self.product] = feature.qualifiers.get('product', [unknown_string])[0]
        return result

def find_frequency_files(directories: List[str]) -> List[str]:
    """Finds all freqs.tsv and masked.tsv files in a given directory."""
    files = []
    for directory in directories:
        for root, _, filenames in os.walk(directory):
            for filename in filenames:
                if filename.endswith((FREQS_FILE_SUFFIX, IVAR_FILE_SUFFIX)):
                    files.append(os.path.join(root, filename))
    return files

def convert_codon_dict_to_string(position_to_nucleotide_dict: Dict[int, str]) -> str:
    """Converts a codon dictionary to a string."""
    return ''.join(position_to_nucleotide_dict[position] for position in sorted(position_to_nucleotide_dict))

def get_reference_codon(codon_positions: Tuple[int, int, int], freqs_data: pd.DataFrame) -> Dict[int, str]:
    """Extracts the reference codon from freqs data based on the position trio."""
    return {position: freqs_data.loc[freqs_data['ref_pos'] == position, 'ref_base'].values[0] for position in codon_positions}