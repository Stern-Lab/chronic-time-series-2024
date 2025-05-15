import re
import requests
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Set
from analyses.mutation_filtering.core.mutation import Mutation

def extract_position(mutation_str: str) -> int:
    """Extract numeric position from a mutation string like 'A123T'."""
    match = re.search(r'\d+', str(mutation_str))
    return int(match.group()) if match else None

def load_problematic_sites(url: str = "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf") -> Set[int]:
    """Download and parse problematic mutation positions from VCF."""
    response = requests.get(url)
    response.raise_for_status()
    lines = [line for line in response.text.splitlines() if not line.startswith("#")]
    return {int(line.split("\t")[1]) for line in lines}

def load_suspicious_mutations(csv_path: Path) -> pd.DataFrame:
    """Load suspicious mutations table and extract relevant fields."""
    df = pd.read_csv(csv_path)
    df["POS"] = df["mut_nuc"].apply(extract_position)
    df["mutation"] = df["mut_nuc"]
    return df

def apply_soft_filter_df_with_reasons(
    unfiltered_mutations: List[Mutation],
    suspicious_set: Set[str],
    problematic_set: Set[int]
) -> Tuple[List[Mutation], pd.DataFrame]:
    """
    Soft filtering with multi-reason support (vectorized).
    Returns:
        - passed: List of Mutation
        - failed_df: DataFrame with semicolon-separated 'drop_reason'
    """
    df = pd.DataFrame([m.to_dict() for m in unfiltered_mutations])

    # Create masks for each filtering reason
    mask_N_ref = df["REF"].str.contains("N", case=False, na=False)
    mask_low_cov = (df["ALT_DP"] < 50) | (df["TOTAL_DP"] < 100)
    mask_suspicious = df["mutation"].isin(suspicious_set)
    mask_problematic = df["POS"].isin(problematic_set)

    # Collect reasons into a list of lists
    reasons_matrix = []
    for i in range(len(df)):
        reasons = []
        if mask_N_ref.iat[i]: reasons.append("N_ref")
        if mask_low_cov.iat[i]: reasons.append("low_coverage")
        if mask_suspicious.iat[i]: reasons.append("suspicious")
        if mask_problematic.iat[i]: reasons.append("problematic")
        reasons_matrix.append(";".join(reasons))

    df["drop_reason"] = reasons_matrix

    # Separate passed and failed
    passed_df = df[df["drop_reason"] == ""].copy()
    failed_df = df[df["drop_reason"] != ""].copy()

    passed = [Mutation(**row) for row in passed_df.to_dict(orient="records")]
    return passed, failed_df


