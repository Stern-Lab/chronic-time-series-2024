# code by Noam
# positions_problematic_virological
# https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473

import io
import pandas as pd


def read_vcf(path=r"./Filter_Usecase/filter_replicate_script/problematic_sites_sarsCov2.vcf"):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': float, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t').rename(columns={'#CHROM': 'CHROM'})


def list_PROBLEMATIC_positions():
    problematic_positions = read_vcf()
    CAUTION_POSITIONS = problematic_positions[problematic_positions.FILTER == 'caution'].POS.tolist()
    MASK_POSITIONS = problematic_positions[problematic_positions.FILTER == 'mask'].POS.tolist()
    HIGHLY_HOMOPLASIC = [187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700]
    HYPERMUTABILIITY = [11083, 15324, 21575]
    PROBLEMATIC = CAUTION_POSITIONS + MASK_POSITIONS + HIGHLY_HOMOPLASIC + HYPERMUTABILIITY
    return PROBLEMATIC