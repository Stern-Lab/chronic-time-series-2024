from dataclasses import dataclass, asdict
from typing import Optional
import pandas as pd

@dataclass
class Mutation:
    POS: int
    mutation: str
    ALT_DP: int
    TOTAL_DP: int
    ALT_FREQ: float
    ALT: str
    REF: str
    ALT_QUAL: float
    REF_QUAL: Optional[float]
    REF_CODON: Optional[str]
    REF_AA: Optional[str]
    mutation_type: Optional[str]
    protein: Optional[str]
    drop_reason: Optional[str] = None

    @classmethod
    def from_row(cls, row: pd.Series) -> "Mutation":
        def is_valid(value):
            return value not in [None, "", "NA"] and pd.notna(value)

        def safe_get(key: str, required: bool = False):
            value = row.get(key)
            if is_valid(value):
                return value
            if required:
                raise ValueError(f"Missing or invalid required field: {key}")
            return None

        pos = int(safe_get("POS", True))
        ref = str(safe_get("REF", True))
        alt = str(safe_get("ALT", True))
        mutation_str = f"{ref}{pos}{alt}"

        ref_aa = safe_get("REF_AA")
        alt_aa = safe_get("ALT_AA")
        mutation_type = f"{ref_aa}{pos}{alt_aa}" if all(is_valid(x) for x in [ref_aa, alt_aa]) else None

        return cls(
            POS=pos,
            mutation=mutation_str,
            ALT_DP=int(safe_get("ALT_DP", True)),
            TOTAL_DP=int(safe_get("TOTAL_DP", True)),
            ALT_FREQ=float(safe_get("ALT_FREQ", True)),
            ALT=alt,
            REF=ref,
            ALT_QUAL=float(safe_get("ALT_QUAL", True)),
            REF_QUAL=float(safe_get("REF_QUAL")) if is_valid(safe_get("REF_QUAL")) else None,
            REF_CODON=safe_get("REF_CODON"),
            REF_AA=ref_aa,
            mutation_type=mutation_type,
            protein=safe_get("GFF_FEATURE")
        )

    def to_dict(self) -> dict:
        return asdict(self)
