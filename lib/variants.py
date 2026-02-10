"""
Variant handling utilities: standardization, ID generation, preparation.
"""

import pandas as pd

from .log import log


def make_pos_key(df: pd.DataFrame, chrom_col: str = 'chrom', pos_col: str = 'pos') -> pd.Series:
    """Create position keys (chrom_pos) for matching."""
    return df[chrom_col].astype(str) + '_' + df[pos_col].astype(str)


def make_variant_id(
    df: pd.DataFrame,
    chrom_col: str = 'chrom',
    pos_col: str = 'pos',
    ref_col: str = 'ref',
    alt_col: str = 'alt',
) -> pd.Series:
    """Create standardized variant IDs (chrom:pos_ref>alt)."""
    return (
        df[chrom_col].astype(str) + ':' +
        df[pos_col].astype(str) + '_' +
        df[ref_col].astype(str) + '>' +
        df[alt_col].astype(str)
    )


def prepare_variants(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize a variants DataFrame for downstream processing.

    - Drops rows with missing essential columns (chrom, pos, ref, alt)
    - Converts pos to int
    - Ensures chr prefix on chromosome names
    - Creates variant_id if not present
    """
    df = df.copy()

    required_cols = ['chrom', 'pos', 'ref', 'alt']
    n_before = len(df)
    df = df.dropna(subset=required_cols)
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        log(f"  Dropped {n_dropped} variants with missing coordinates")

    df['pos'] = df['pos'].astype(int)

    if len(df) > 0 and not str(df['chrom'].iloc[0]).startswith('chr'):
        df['chrom'] = 'chr' + df['chrom'].astype(str)

    if 'variant_id' not in df.columns:
        df['variant_id'] = make_variant_id(df)

    return df
