import io
import pandas as pd


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': float, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t').rename(columns={'#CHROM': 'CHROM'})


def list_problematic_sites():
    problematic_sites = read_vcf('./input/problematic_sites_sarsCov2.vcf')
    CAUTION_POSITIONS = problematic_sites[problematic_sites.FILTER == 'caution'].POS.tolist()
    MASK_POSITIONS = problematic_sites[problematic_sites.FILTER == 'mask'].POS.tolist()
    HIGHLY_HOMOPLASIC = [187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700]
    HYPERMUTABILIITY = [11083, 15324, 21575]
    PROBLEMATIC = CAUTION_POSITIONS + MASK_POSITIONS + HIGHLY_HOMOPLASIC + HYPERMUTABILIITY
    return PROBLEMATIC


def phase1_filtering(coverage, frequency, base_count, coverage_t=100, frequency_t=0.01, base_count_t=50):
    """"""
    if coverage < coverage_t:
        return None

    else:
        if frequency < frequency_t:
            return 0

        else:
            if base_count < base_count_t:
                return None

            else:
                return frequency


def calc_freq_weighted_avg(base_cnt1, base_cnt2, cov1, cov2):
    weighted_avg = (base_cnt1 + base_cnt2) / (cov1 + cov2)
    return weighted_avg


def phase2_filtering_technical_replicates(freq1, freq2, base_cnt1, base_cnt2, cov1, cov2):
    """"""
    if pd.isna(freq1) or pd.isna(freq2):
        return None

    else:
        if freq1 == 0 and freq2 == 0:
            return 0

        elif (freq1 == 0 and freq2 != 0) or (freq1 != 0 and freq2 == 0):
            return None

        else:  # f1!=0 and f2!=0
            if freq1 < 0.5 and freq2 < 0.5:
                delta_lim = 0.1
            else:
                delta_lim = 0.3

            delta = abs(freq1 - freq2)

            if delta > delta_lim:
                return None
            else:
                return calc_freq_weighted_avg(base_cnt1, base_cnt2, cov1, cov2)


def filter_mutations(df, coverage_t: int, frequency_t: float, base_count_t: int, ignore_indels: bool):
    """
    Gets a freq Dataframe, and filter it by the desired coverage, frequency and base count thresholds (in that order).
    :return: a filtered DataFrame
    """

    # filter out problematic positions
    problematic = list_problematic_sites()
    filtered_df = df[~df['POS'].astype(int).isin(problematic)].copy()

    # phase 0 - get only mutations - unnecessary since iVar outputs only variants

    # filter
    if ignore_indels:
        filtered_df = filtered_df.loc[(filtered_df['ALT'].str.len() == 1)]  # indels have + or - before the bases

    filtered_df['new_ALT_FREQ'] = filtered_df.apply(lambda row: phase1_filtering(row['TOTAL_DP'], row['ALT_FREQ'],
                                                                                 row['ALT_DP'],
                                                                                 coverage_t=coverage_t,
                                                                                 frequency_t=frequency_t,
                                                                                 base_count_t=base_count_t), axis=1)

    filtered_df.dropna(subset=['new_ALT_FREQ'], inplace=True)

    return filtered_df


def get_technical_reps(df, y_df=None):
    """"""
    if y_df:    # allows to manually work with two DataFrames
        x_df = df

    else:
        # split into replicates
        x_df = df.loc[(df['replicate'] == 'A')]
        y_df = df.loc[(df['replicate'] == 'B')]     # fixed from y_df = df.loc[(df['replicate'] != 'B')]

    # Merge DataFrame based on sample and mutation
    merged_df = pd.merge(x_df[['sample', 'mutation', 'POS', 'ALT', 'REF', 'transition', 'type', 'source',
                               'ALT_DP', 'REF_QUAL', 'ALT_QUAL', 'TOTAL_DP', 'PVAL', 'new_ALT_FREQ',
                               'PVAL', 'PASS', 'GFF_FEATURE', 'REF_CODON', 'REF_AA', 'ALT_CODON', 'ALT_AA']],
                         y_df[['sample', 'mutation', 'source',
                               'ALT_DP', 'REF_QUAL', 'ALT_QUAL', 'TOTAL_DP', 'PVAL', 'PASS', 'new_ALT_FREQ']],
                         left_on=['sample', 'mutation'], right_on=['sample', 'mutation'], how='inner')

    # get new frequencies based on filtering decision tree phase 2
    merged_df['merged_ALT_DP'] = merged_df['ALT_DP_x'] + merged_df['ALT_DP_y']
    merged_df['merged_ALT_QUAL'] = ((merged_df['ALT_QUAL_x']*merged_df['ALT_DP_x']) + (merged_df['ALT_QUAL_y']*merged_df['ALT_DP_y'])) / (merged_df['ALT_DP_x'] + merged_df['ALT_DP_y'])
    merged_df['merged_TOTAL_DP'] = merged_df['TOTAL_DP_x'] + merged_df['TOTAL_DP_y']
    merged_df['ALT_FREQ_delta'] = abs(merged_df['new_ALT_FREQ_x'] - merged_df['new_ALT_FREQ_y'])
    merged_df['merged_ALT_FREQ'] = merged_df.apply(lambda row:
                                                    phase2_filtering_technical_replicates(row['new_ALT_FREQ_x'],
                                                                                          row['new_ALT_FREQ_y'],
                                                                                          row['ALT_DP_x'],
                                                                                          row['ALT_DP_y'],
                                                                                          row['TOTAL_DP_x'],
                                                                                          row['TOTAL_DP_y']), axis=1)
    # drop freqs that were deemed NA
    filtered_merged_df = merged_df.dropna(subset=['merged_ALT_FREQ'])

    return filtered_merged_df

