import pandas as pd
from Bio import SeqIO
import glob


def assign_transition_type(t):
    transitions = ['AG', 'GA', 'TC', 'CT']
    oxidation = ['CA', 'GT']
    transversions = ['AC', 'TG', 'TA', 'AT', 'GC', 'CG']
    if not isinstance(t, str):
        return 'err'
    if t in transitions:
        return 'ts'
    elif t in oxidation:
        return 'ox'
    elif t[0] == '-':
        return 'ins'
    elif t[1] == '-':
        return 'del'
    elif t in transversions:
        return 'tv'
    else:
        return 'ref'


def add_data_columns(df):
    df['sample'] = df.source.map(lambda s: s.split('-')[0])     # get the sample id
    df['replicate'] = df.source.map(lambda s: s.split('-')[1])  # get replicate id
    df['mutation'] = df['REF'] + df['POS'].astype(int).map(str) + df['ALT']
    df['transition'] = df['REF'] + df['ALT']
    df['type'] = df.transition.map(assign_transition_type)


def merge_variants_w_info(var_df, sym_onset_path, ct_values_path):
    # read info tables
    sym_onset = pd.read_csv(sym_onset_path)
    ct_values = pd.read_csv(ct_values_path)

    # merge
    merged_w_sym_onset = pd.merge(var_df, sym_onset, how='left', on=['sample'])
    merged_variants = pd.merge(merged_w_sym_onset, ct_values[["sample", "N_gene_ct"]], how='left', on=['sample'])

    # drop unnecessary columns & return df
    merged_variants.dropna(how='all', axis='columns')
    return merged_variants


def read_gff_file(gff_file_path: str):
    """
    Read a GFF file and extract CDS regions.

    Parameters:
    gff_file (str): Path to the GFF file.

    Returns:
    dict: A dictionary with gene IDs as keys and list of CDS regions as values.
    """
    cds_regions = {}
    with open(gff_file_path) as in_handle:
        for line in in_handle:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if parts[2] == "CDS":
                    seqname, source, feature, start, end, score, strand, frame, attributes = parts
                    start, end = int(start), int(end)
                    strand = 1 if strand == '+' else -1
                    gene_id = attributes.split('=')[1]
                    if gene_id not in cds_regions:
                        cds_regions[gene_id] = []
                    cds_regions[gene_id].append((start, end, strand))
    return cds_regions


def read_bed_file(bed_file):
    """
    Reads a BED file and returns a DataFrame with the primer information.

    Parameters:
    bed_file (str): Path to the BED file.

    Returns:
    pd.DataFrame: DataFrame containing the primer information.
    """
    columns = ["chrom", "start", "end", "name", "score", "strand", "sequence"]
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=columns)
    return bed_df


def nucleotide_to_aa_position(nucleotide_position: int, gff_file_path: str, gene_id: str = None):
    """
    Convert nucleotide position to amino acid position using GFF file.

    Parameters:
    nucleotide_position (int): The position of the nucleotide.
    gff_file (str): Path to the GFF file.
    gene_id (str): The gene ID to consider (optional). If None, considers all genes.

    Returns:
    int, str: The amino acid position and the relevant gene ID.
    """
    # Read CDS regions from the GFF file
    cds_regions = read_gff_file(gff_file_path)

    if gene_id and gene_id in cds_regions:
        gene_cds_regions = {gene_id: cds_regions[gene_id][0]}
    else:
        # If no gene_id provided, search through all genes
        gene_cds_regions = {gene: cds_list[0] for gene, cds_list in cds_regions.items()}

    # Find the relevant CDS region
    cumulative_length = 0
    for gene, (start, end, strand) in gene_cds_regions.items():
        start, end, nucleotide_position = int(start), int(end), int(nucleotide_position)
        if start <= nucleotide_position <= end:
            relative_position = cumulative_length + (nucleotide_position - start + 1)
            aa_position = (relative_position - 1) // 3 + 1
            return aa_position, gene
        else:
            cumulative_length += (end - start + 1)

    return None, None     # if position is not in CDS


def categorize_synonymity(row):
    if row['REF_AA'] == row['ALT_AA']:
        return 'syn'
    elif row['ALT_AA'] == '*':
        return 'stop'
    else:
        return 'non-syn'


def categorize_cds(row, gff_file_path: str):
    if row['GFF_FEATURE']:
        return 'CDS'
    else:
        cds_regions = read_gff_file(gff_file_path)
        first_gene = list(cds_regions.keys())[0]
        first_gene_start = cds_regions[first_gene][0][0]
        last_gene = list(cds_regions.keys())[-1]
        last_gene_end = cds_regions[last_gene][-1][1]
        if row['POS'] < first_gene_start:
            return '5UTR'
        elif row['POS'] > last_gene_end:
            return '3UTR'
        else:
            return 'TRS'


def get_GC_content_near_pos(pos: int, consensus_seq_path: str, window_size: int = 25):
    """
    Calculate the GC content in a window around a given position.

    Parameters:
    pos (int): The position of the nucleotide.
    cons_seq (str): The consensus sequence.
    window_size (int): The size of the window around the position.

    Returns:
    float: The GC content in the window around the position.
    """

    # Read the consensus sequence
    with open(consensus_seq_path) as in_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            cons_seq = str(record.seq)

    # Get the start and end positions of the window
    start_pos = max(0, pos - window_size)
    end_pos = min(len(cons_seq), pos + window_size)

    # Get the sequence in the window
    window_seq = cons_seq[int(start_pos):int(end_pos)]

    # Calculate the GC content
    gc_content = (window_seq.count('G') + window_seq.count('C')) / len(window_seq)
    gc_perc = gc_content * 100
    return gc_perc


def get_distance_from_closest_primer(nt_pos: int, primers_version: str, primer_data_dir_path: str):
    """
    Calculates the distance from the nucleotide position to the nearest primer.

    Parameters:
    nt_pos (int): The position of the nucleotide.
    primers_version (str): The version of the primers. Should be one of 'V1', 'V2', 'V3', 'V4', 'V4.1', or 'V5.3.2'.
    primer_data_dir_path (str): Path to the directory containing the BED files.

    Returns:
    tuple: The distance to the nearest primer (0 if inside the primer's range) and the strand ('+' for forward, '-' for reverse).
    """
    if primers_version.upper() not in ['V1', 'V2', 'V3', 'V4', 'V4.1', 'V5.3.2']:
        raise ValueError("Invalid primers version. Please choose from V1, V2, V3, V4, V4.1, or V5.3.2")
    else:
        # get the path to the primer BED file & read it
        bed_file = glob.glob(f'{primer_data_dir_path}/{primers_version}/*.primer.bed')  # this is a list of 1 str
        bed_df = read_bed_file(bed_file[0])

        # convert the nucleotide position from str to int
        nt_pos = int(nt_pos)

        # Calculate the distances
        min_distance = float('inf')
        closest_strand = None

        for index, row in bed_df.iterrows():
            if row['start'] <= nt_pos <= row['end']:
                # If the nucleotide position is within the primer's range, distance is 0
                return 0, row['strand']
            else:
                # Calculate the distance to the start and end of the primer
                distance_to_start = abs(nt_pos - row['start'])
                distance_to_end = abs(nt_pos - row['end'])
                min_dist = min(distance_to_start, distance_to_end)

                if min_dist < min_distance:
                    min_distance = min_dist
                    closest_strand = row['strand']

        # Return the minimum distance and the strand
        return min_distance, closest_strand


def aa_to_orf1ab_position(aa_position: int, region: str):
    """
    Converts an amino acid position within the ORF1a or ORF1b region to its position within the ORF1ab region.

    Parameters:
    aa_position (int): The amino acid position within the ORF1a or ORF1b region.
    region (str): The region name ('ORF1a' or 'ORF1b').

    Returns:
    int: The amino acid position within the ORF1ab region.
    """
    # Lengths of ORF1a and ORF1b in amino acids
    length_orf1a = 4401  # Length of ORF1a in amino acids
    # ORF1b follows ORF1a, so the start position of ORF1b in ORF1ab is length_orf1a + 1

    if region == 'ORF1a':
        return aa_position
    elif region == 'ORF1b':
        return length_orf1a + aa_position
    else:
        raise ValueError("Invalid region name. Use 'ORF1a' or 'ORF1b'.")


def adjust_orf1ab_pos_aa_and_gff_feature(df: pd.DataFrame):
    """
    Adjusts the POS_AA and GFF_FEATURE columns in the DataFrame for rows where GFF_FEATURE is either 'ORF1a' or 'ORF1b'.

    Parameters:
    df (pd.DataFrame): The input DataFrame.

    Returns:
    pd.DataFrame: The adjusted DataFrame.
    """
    def adjust_row(row):
        if row['GFF_FEATURE'] in ['ORF1a', 'ORF1b']:
            row['ORF1'] = row['GFF_FEATURE']
            row['POS_AA'] = aa_to_orf1ab_position(aa_position=row['POS_AA'], region=row['GFF_FEATURE'])
            row['GFF_FEATURE'] = 'ORF1ab'
        else:
            row['ORF1'] = None
        return row

    # Apply the adjustment to each row
    df = df.apply(adjust_row, axis=1)
    return df
