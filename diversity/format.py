COL_REF_POS = "ref_pos"
COL_READ_BASE = "read_base"
COL_BASE_COUNT = "base_count"
COL_COVERAGE = "coverage"
COL_FREQUENCY = "frequency"
COL_BASE_RANK = "base_rank"


class FrequencyTableRow:
    REQUIRED_COLUMNS = [
        COL_REF_POS,
        COL_READ_BASE,
        COL_BASE_COUNT,
        COL_COVERAGE,
        COL_FREQUENCY,
        COL_BASE_RANK
    ]

    @staticmethod
    def verify_columns(df):
        missing = [col for col in FrequencyTableRow.REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")