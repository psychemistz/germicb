"""
Variant matching against eQTL reference databases and concordance computation.
"""

from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy import stats

from .log import log
from .variants import make_pos_key


def match_by_rsid(
    variants_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    ref_name: str,
    beta_col: str = 'beta',
    pval_col: str = 'pval',
    extra_cols: Optional[list] = None,
) -> pd.DataFrame:
    """
    Match variants to a reference eQTL database by rsID.

    Use this instead of match_by_position() when variant coordinates are on
    different genome builds (e.g. therapy variants on hg38, OneK1K on hg19).

    Args:
        variants_df: Therapy variants DataFrame (must have rsid column)
        ref_df: Reference eQTL DataFrame (must have rsid column)
        ref_name: Short name for the reference (e.g. 'onek1k')
        beta_col: Column name for effect size in ref_df
        pval_col: Column name for p-value in ref_df
        extra_cols: Additional columns to carry from ref_df

    Returns:
        variants_df with added {ref_name}_* columns and {ref_name}_matched boolean
    """
    log(f"Matching to {ref_name} eQTLs (by rsID)...")

    variants_df = variants_df.copy()
    ref_df = ref_df.copy()

    # Drop rows without valid rsIDs
    ref_valid = ref_df[ref_df['rsid'].notna() & (ref_df['rsid'] != '') & (ref_df['rsid'] != '.')].copy()
    var_valid_mask = variants_df['rsid'].notna() & (variants_df['rsid'] != '') & (variants_df['rsid'] != '.')

    common_rsids = set(variants_df.loc[var_valid_mask, 'rsid']) & set(ref_valid['rsid'])
    log(f"  Matching rsIDs: {len(common_rsids):,}")

    if len(common_rsids) == 0:
        variants_df[f'{ref_name}_matched'] = False
        return variants_df

    # Best hit per rsID (lowest p-value)
    ref_best = ref_valid.loc[ref_valid.groupby('rsid')[pval_col].idxmin()]

    # Build merge columns
    merge_cols = ['rsid', beta_col, pval_col]
    rename_map = {
        beta_col: f'{ref_name}_beta',
        pval_col: f'{ref_name}_pval',
    }

    # Detect gene column
    for gene_col in ('gene_symbol', 'gene'):
        if gene_col in ref_best.columns:
            merge_cols.append(gene_col)
            rename_map[gene_col] = f'{ref_name}_gene'
            break

    # Detect celltype column
    for ct_col in (f'{ref_name}_celltype', 'celltype', 'cell_type'):
        if ct_col in ref_best.columns:
            merge_cols.append(ct_col)
            if ct_col != f'{ref_name}_celltype':
                rename_map[ct_col] = f'{ref_name}_celltype'
            break

    if extra_cols:
        for col in extra_cols:
            if col in ref_best.columns and col not in merge_cols:
                merge_cols.append(col)

    ref_subset = ref_best[merge_cols].copy()
    ref_subset = ref_subset.rename(columns=rename_map)

    variants_df = variants_df.merge(ref_subset, on='rsid', how='left')
    variants_df[f'{ref_name}_matched'] = variants_df[f'{ref_name}_beta'].notna()

    n_matched = variants_df[f'{ref_name}_matched'].sum()
    log(f"  Matched: {n_matched:,} / {len(variants_df):,} ({100*n_matched/len(variants_df):.1f}%)")

    return variants_df


def match_by_position(
    variants_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    ref_name: str,
    beta_col: str = 'beta',
    pval_col: str = 'pval',
    extra_cols: Optional[list] = None,
) -> pd.DataFrame:
    """
    Match variants to a reference eQTL database by genomic position.

    Args:
        variants_df: Therapy variants DataFrame (must have chrom, pos columns)
        ref_df: Reference eQTL DataFrame (must have chrom, pos columns)
        ref_name: Short name for the reference (e.g. 'dice', 'onek1k')
        beta_col: Column name for effect size in ref_df
        pval_col: Column name for p-value in ref_df
        extra_cols: Additional columns to carry from ref_df

    Returns:
        variants_df with added {ref_name}_* columns and {ref_name}_matched boolean
    """
    log(f"Matching to {ref_name} eQTLs...")

    variants_df = variants_df.copy()
    variants_df['pos_key'] = make_pos_key(variants_df)

    ref_df = ref_df.copy()
    ref_df['pos_key'] = make_pos_key(ref_df)

    common_positions = set(variants_df['pos_key']) & set(ref_df['pos_key'])
    log(f"  Matching positions: {len(common_positions)}")

    if len(common_positions) == 0:
        variants_df[f'{ref_name}_matched'] = False
        return variants_df

    # Best hit per position (lowest p-value)
    ref_best = ref_df.loc[ref_df.groupby('pos_key')[pval_col].idxmin()]

    # Build merge columns
    merge_cols = ['pos_key', beta_col, pval_col]
    rename_map = {
        beta_col: f'{ref_name}_beta',
        pval_col: f'{ref_name}_pval',
    }

    # Detect gene column
    for gene_col in ('gene_symbol', 'gene'):
        if gene_col in ref_best.columns:
            merge_cols.append(gene_col)
            rename_map[gene_col] = f'{ref_name}_gene'
            break

    # Detect celltype column
    for ct_col in (f'{ref_name}_celltype', 'celltype', 'cell_type'):
        if ct_col in ref_best.columns:
            merge_cols.append(ct_col)
            if ct_col != f'{ref_name}_celltype':
                rename_map[ct_col] = f'{ref_name}_celltype'
            break

    if extra_cols:
        for col in extra_cols:
            if col in ref_best.columns and col not in merge_cols:
                merge_cols.append(col)

    ref_subset = ref_best[merge_cols].copy()
    ref_subset = ref_subset.rename(columns=rename_map)

    variants_df = variants_df.merge(ref_subset, on='pos_key', how='left')
    variants_df[f'{ref_name}_matched'] = variants_df[f'{ref_name}_beta'].notna()

    n_matched = variants_df[f'{ref_name}_matched'].sum()
    log(f"  Matched: {n_matched:,} / {len(variants_df):,} ({100*n_matched/len(variants_df):.1f}%)")

    return variants_df


def compute_concordance(df: pd.DataFrame, col_a: str, col_b: str) -> Dict:
    """
    Compute direction concordance and Spearman correlation between two columns.

    Args:
        df: DataFrame containing both columns
        col_a: First column name (e.g. AlphaGenome impact)
        col_b: Second column name (e.g. eQTL beta)

    Returns:
        Dict with n_valid, direction_concordance, spearman_r, spearman_pval
    """
    valid = df[
        df[col_a].notna() & df[col_b].notna() &
        np.isfinite(df[col_a]) & np.isfinite(df[col_b]) &
        (df[col_a] != 0) & (df[col_b] != 0)
    ]

    if len(valid) == 0:
        return {
            'n_valid': 0,
            'direction_concordance': 0.0,
            'spearman_r': 0.0,
            'spearman_pval': 1.0,
        }

    concordance = float((np.sign(valid[col_a]) == np.sign(valid[col_b])).mean())
    corr, pval = stats.spearmanr(valid[col_a].values, valid[col_b].values, nan_policy='omit')

    return {
        'n_valid': len(valid),
        'direction_concordance': concordance,
        'spearman_r': float(corr) if not np.isnan(corr) else 0.0,
        'spearman_pval': float(pval) if not np.isnan(pval) else 1.0,
    }
