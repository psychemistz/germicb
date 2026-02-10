#!/usr/bin/env python3
"""
Validate AlphaGenome Predictions Against eQTL Databases
========================================================
Prioritizing Immunotherapy-Associated Regulatory Germline Variants Through AlphaGenome

Match immunotherapy-associated variants against independent eQTL sources:

1. DICE - Immune cell type-specific eQTLs (12 cell types, hg38)
2. OneK1K - Single-cell PBMC eQTLs (14 cell types, hg38; Yazar et al., Science 2022)
3. Cohort eQTLs - Per-cohort local eQTL analysis (sanity check, hg19→hg38 liftover)

GTEx is EXCLUDED from validation — AlphaGenome is trained on GTEx RNA-seq data,
making GTEx concordance circular (memorized training signal, not predictive power).

Input:
- results/alphagenome/alphagenome_predictions.h5ad
- DICE hg38 data (from cytokine atlas liftover)
- OneK1K sc-eQTLs (data/eqtl_references/oneK1K/esnp_table.tsv.gz)
- Per-cohort 9eQTL/output/fin_eqtl_tbl.txt (optional sanity check)

Output:
- results/eqtl_validation/eqtl_matched.csv
- results/eqtl_validation/validation_metrics.json
- results/eqtl_validation/validation_report.md
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import anndata as ad

# ==============================================================================
# Configuration
# ==============================================================================

PROJECT_DIR = Path('/data/parks34/projects/4germicb')
DATA_DIR = PROJECT_DIR / 'data'
RESULTS_DIR = PROJECT_DIR / 'results'
ALPHAGENOME_DIR = RESULTS_DIR / 'alphagenome'
OUTPUT_DIR = RESULTS_DIR / 'eqtl_validation'

# DICE hg38 (lifted from hg19)
DICE_HG38_DIR = PROJECT_DIR / 'data' / 'eqtl_references' / 'dice' / 'hg38'

# OneK1K sc-eQTLs (hg38, Yazar et al. Science 2022)
ONEK1K_PATH = PROJECT_DIR / 'data' / 'eqtl_references' / 'oneK1K' / 'esnp_table.tsv.gz'

# CIMA eQTLs (optional)
CIMA_EQTL_DIR = Path('/data/Jiang_Lab/Data/Seongyong/CIMA/xQTL')

# Cohort eQTL directories (for sanity check)
WES_DIR = Path('/data/parks34/projects/project_outdated/WES_ImmunoPredict/FINISHED')


def log(msg: str):
    """Print timestamped log message."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_alphagenome_predictions(input_path: Path) -> pd.DataFrame:
    """
    Load AlphaGenome predictions from h5ad file.

    Returns DataFrame with variant info and track predictions.
    """
    log(f"Loading AlphaGenome predictions: {input_path}")

    adata = ad.read_h5ad(input_path)
    log(f"  Loaded {adata.n_obs:,} variants, {adata.n_vars} tracks")

    # Get variant metadata
    variants_df = adata.obs.copy()

    # Compute summary impact scores from prediction matrix
    X = adata.X

    # Mean absolute impact across all tracks
    variants_df['alphagenome_impact_mean'] = np.abs(X).mean(axis=1)
    variants_df['alphagenome_impact_max'] = np.abs(X).max(axis=1)

    # Track-specific impacts (if available)
    track_names = adata.var_names.tolist()
    for i, track in enumerate(track_names):
        clean_name = track.replace('RNA_SEQ_GTEx_', '').replace(' ', '_')
        variants_df[f'ag_{clean_name}'] = X[:, i]

    log(f"  Impact score range: {variants_df['alphagenome_impact_mean'].min():.4f} to {variants_df['alphagenome_impact_mean'].max():.4f}")

    return variants_df


def load_dice_eqtls() -> Optional[pd.DataFrame]:
    """
    Load DICE immune cell eQTLs (hg38 lifted coordinates).

    Returns DataFrame with columns: chrom, pos, ref, alt, gene_symbol, beta, pval, dice_celltype
    """
    if not DICE_HG38_DIR.exists():
        log(f"  DICE hg38 directory not found: {DICE_HG38_DIR}")
        return None

    log(f"Loading DICE eQTLs from: {DICE_HG38_DIR}")

    files_to_load = list(DICE_HG38_DIR.glob('*_hg38.tsv'))
    if not files_to_load:
        log("  No DICE hg38 TSV files found")
        return None

    log(f"  Found {len(files_to_load)} cell type files")

    all_records = []

    for file_path in files_to_load:
        if file_path.stat().st_size == 0:
            continue

        dice_celltype = file_path.stem.replace('_hg38', '')

        df = pd.read_csv(file_path, sep='\t')

        for _, row in df.iterrows():
            # Parse INFO field
            info_dict = {}
            if pd.notna(row.get('INFO', '')):
                for item in str(row['INFO']).split(';'):
                    if '=' in item:
                        key, val = item.split('=', 1)
                        info_dict[key] = val

            # Use hg38 coordinates
            if pd.notna(row.get('CHROM_hg38')) and pd.notna(row.get('POS_hg38')):
                all_records.append({
                    'chrom': row['CHROM_hg38'],
                    'pos': int(row['POS_hg38']),
                    'rsid': row.get('ID', ''),
                    'ref': row.get('REF', ''),
                    'alt': row.get('ALT', ''),
                    'gene_symbol': info_dict.get('GeneSymbol', ''),
                    'gene_id': info_dict.get('Gene', ''),
                    'beta': float(info_dict.get('Beta', 0)),
                    'pval': float(info_dict.get('Pvalue', 1)),
                    'dice_celltype': dice_celltype,
                })

    if not all_records:
        log("  No DICE records loaded")
        return None

    dice_df = pd.DataFrame(all_records)
    log(f"  Total DICE eQTLs: {len(dice_df):,} from {dice_df['dice_celltype'].nunique()} cell types")

    return dice_df


def load_onek1k_eqtls() -> Optional[pd.DataFrame]:
    """
    Load OneK1K single-cell PBMC eQTLs (hg38).

    Source: Yazar et al., Science 2022 — ~1.27M PBMCs from 982 donors, 14 cell types.
    Independent of AlphaGenome training data.

    Returns DataFrame with columns: chrom, pos, rsid, gene, beta, pval, onek1k_celltype
    """
    if not ONEK1K_PATH.exists():
        log(f"  OneK1K file not found: {ONEK1K_PATH}")
        return None

    log(f"Loading OneK1K eQTLs: {ONEK1K_PATH}")

    onek1k_df = pd.read_csv(ONEK1K_PATH, sep='\t')
    log(f"  Loaded {len(onek1k_df):,} records")

    # Filter to significant eQTLs (FDR < 0.05)
    if 'FDR' in onek1k_df.columns:
        onek1k_df = onek1k_df[onek1k_df['FDR'] < 0.05].copy()
        log(f"  After FDR < 0.05 filter: {len(onek1k_df):,}")

    # Standardize columns
    # CHR column has no "chr" prefix in OneK1K
    onek1k_df['chrom'] = 'chr' + onek1k_df['CHR'].astype(str)
    onek1k_df['pos'] = onek1k_df['POS'].astype(int)
    onek1k_df['rsid'] = onek1k_df['RSID']
    onek1k_df['gene'] = onek1k_df['GENE']
    onek1k_df['beta'] = onek1k_df['SPEARMANS_RHO'].astype(float)
    onek1k_df['pval'] = onek1k_df['P_VALUE'].astype(float)
    onek1k_df['onek1k_celltype'] = onek1k_df['CELL_TYPE']

    log(f"  Cell types: {onek1k_df['onek1k_celltype'].nunique()}")
    log(f"  Unique variants: {onek1k_df['rsid'].nunique():,}")

    return onek1k_df[['chrom', 'pos', 'rsid', 'gene', 'beta', 'pval', 'onek1k_celltype']].copy()


def load_cohort_eqtls() -> Optional[pd.DataFrame]:
    """
    Load per-cohort local eQTL results for sanity check.

    These are from the WES immunotherapy cohort analyses (hg19 coordinates).
    Requires liftover to hg38 for matching against DICE/OneK1K/AlphaGenome.

    Returns DataFrame with columns: chrom, pos_hg19, pos_hg38, rsid, gene, beta, pval, fdr, cohort
    """
    if not WES_DIR.exists():
        log(f"  WES directory not found: {WES_DIR}")
        return None

    log(f"Loading cohort eQTLs from: {WES_DIR}")

    # Try importing pyliftover for hg19→hg38 conversion
    try:
        from pyliftover import LiftOver
        lo = LiftOver('hg19', 'hg38')
    except ImportError:
        log("  WARNING: pyliftover not available, using rsid matching only")
        lo = None

    all_cohort_eqtls = []
    cohort_dirs = [d for d in WES_DIR.iterdir()
                   if d.is_dir() and not d.name.startswith(('summary', 'logis', 'ridge', 'Trans'))]

    for cohort_dir in sorted(cohort_dirs):
        eqtl_file = cohort_dir / '9eQTL' / 'output' / 'fin_eqtl_tbl.txt'
        snp_annot_file = cohort_dir / '9eQTL' / 'data' / 'snp_annot.txt'

        if not eqtl_file.exists():
            continue

        # Load eQTL results (space-delimited)
        try:
            eqtl_df = pd.read_csv(eqtl_file, sep=r'\s+')
        except Exception as e:
            log(f"  WARNING: Could not load {eqtl_file}: {e}")
            continue

        if eqtl_df.empty:
            continue

        # Filter to significant eQTLs
        eqtl_df = eqtl_df[eqtl_df['FDR'] < 0.05].copy()
        if eqtl_df.empty:
            continue

        # Clean beta: replace Inf/-Inf with NaN
        eqtl_df['beta'] = pd.to_numeric(eqtl_df['beta'], errors='coerce')
        eqtl_df = eqtl_df[eqtl_df['beta'].notna() & np.isfinite(eqtl_df['beta'])].copy()

        eqtl_df['cohort'] = cohort_dir.name

        # Join with snp_annot.txt to get positions
        if snp_annot_file.exists():
            try:
                snp_annot = pd.read_csv(snp_annot_file, sep='\t')
                # Match by rsid
                snp_annot = snp_annot.rename(columns={'CHROM': 'chrom', 'POS': 'pos_hg19', 'ID': 'rsid_annot'})
                eqtl_df = eqtl_df.merge(
                    snp_annot[['chrom', 'pos_hg19', 'rsid_annot']].drop_duplicates('rsid_annot'),
                    left_on='snps', right_on='rsid_annot', how='inner'
                )
            except Exception as e:
                log(f"  WARNING: Could not join snp_annot for {cohort_dir.name}: {e}")
                continue
        else:
            continue

        if eqtl_df.empty:
            continue

        # Liftover hg19 → hg38
        if lo is not None:
            hg38_positions = []
            for _, row in eqtl_df.iterrows():
                chrom = row['chrom']
                pos = int(row['pos_hg19'])
                try:
                    result = lo.convert_coordinate(chrom, pos)
                    if result and len(result) > 0:
                        hg38_positions.append(int(result[0][1]))
                    else:
                        hg38_positions.append(np.nan)
                except Exception:
                    hg38_positions.append(np.nan)
            eqtl_df['pos_hg38'] = hg38_positions
            eqtl_df = eqtl_df[eqtl_df['pos_hg38'].notna()].copy()
            eqtl_df['pos_hg38'] = eqtl_df['pos_hg38'].astype(int)

        n_eqtls = len(eqtl_df)
        if n_eqtls > 0:
            log(f"  {cohort_dir.name}: {n_eqtls:,} significant eQTLs (FDR<0.05)")
            all_cohort_eqtls.append(eqtl_df)

    if not all_cohort_eqtls:
        log("  No cohort eQTLs loaded")
        return None

    combined = pd.concat(all_cohort_eqtls, ignore_index=True)
    log(f"  Total cohort eQTLs: {len(combined):,} from {combined['cohort'].nunique()} cohorts")

    return combined


def match_to_dice(variants_df: pd.DataFrame, dice_df: pd.DataFrame) -> pd.DataFrame:
    """
    Match therapy variants to DICE eQTLs by position.

    Returns variants_df with added DICE columns.
    """
    log("Matching to DICE eQTLs...")

    # Create position keys
    variants_df = variants_df.copy()
    variants_df['pos_key'] = variants_df['chrom'].astype(str) + '_' + variants_df['pos'].astype(str)

    dice_df = dice_df.copy()
    dice_df['pos_key'] = dice_df['chrom'].astype(str) + '_' + dice_df['pos'].astype(str)

    # Find matches
    common_positions = set(variants_df['pos_key']) & set(dice_df['pos_key'])
    log(f"  Matching positions: {len(common_positions)}")

    if len(common_positions) == 0:
        variants_df['dice_matched'] = False
        return variants_df

    # Get best DICE hit per position (lowest p-value)
    dice_best = dice_df.loc[dice_df.groupby('pos_key')['pval'].idxmin()]

    # Merge
    dice_subset = dice_best[['pos_key', 'gene_symbol', 'beta', 'pval', 'dice_celltype']].copy()
    dice_subset = dice_subset.rename(columns={
        'gene_symbol': 'dice_gene',
        'beta': 'dice_beta',
        'pval': 'dice_pval',
        'dice_celltype': 'dice_celltype',
    })

    variants_df = variants_df.merge(dice_subset, on='pos_key', how='left')
    variants_df['dice_matched'] = variants_df['dice_gene'].notna()

    n_matched = variants_df['dice_matched'].sum()
    log(f"  Matched: {n_matched:,} / {len(variants_df):,} ({100*n_matched/len(variants_df):.1f}%)")

    return variants_df


def match_to_onek1k(variants_df: pd.DataFrame, onek1k_df: pd.DataFrame) -> pd.DataFrame:
    """
    Match therapy variants to OneK1K sc-eQTLs by position.

    Returns variants_df with added OneK1K columns.
    """
    log("Matching to OneK1K eQTLs...")

    # Create position keys
    variants_df = variants_df.copy()
    variants_df['pos_key'] = variants_df['chrom'].astype(str) + '_' + variants_df['pos'].astype(str)

    onek1k_df = onek1k_df.copy()
    onek1k_df['pos_key'] = onek1k_df['chrom'].astype(str) + '_' + onek1k_df['pos'].astype(str)

    # Find matches
    common_positions = set(variants_df['pos_key']) & set(onek1k_df['pos_key'])
    log(f"  Matching positions: {len(common_positions)}")

    if len(common_positions) == 0:
        variants_df['onek1k_matched'] = False
        return variants_df

    # Get best OneK1K hit per position (lowest p-value)
    onek1k_best = onek1k_df.loc[onek1k_df.groupby('pos_key')['pval'].idxmin()]

    # Merge
    onek1k_subset = onek1k_best[['pos_key', 'gene', 'beta', 'pval', 'onek1k_celltype']].copy()
    onek1k_subset = onek1k_subset.rename(columns={
        'gene': 'onek1k_gene',
        'beta': 'onek1k_beta',
        'pval': 'onek1k_pval',
        'onek1k_celltype': 'onek1k_celltype',
    })

    variants_df = variants_df.merge(onek1k_subset, on='pos_key', how='left')
    variants_df['onek1k_matched'] = variants_df['onek1k_gene'].notna()

    n_matched = variants_df['onek1k_matched'].sum()
    log(f"  Matched: {n_matched:,} / {len(variants_df):,} ({100*n_matched/len(variants_df):.1f}%)")

    return variants_df


def run_cohort_sanity_check(
    variants_df: pd.DataFrame,
    cohort_eqtls: pd.DataFrame,
    dice_df: Optional[pd.DataFrame],
    onek1k_df: Optional[pd.DataFrame],
) -> Dict:
    """
    Sanity check: compare per-cohort eQTLs against reference databases and AlphaGenome.

    Three-way concordance:
    1. Cohort eQTL vs DICE — do small-cohort eQTLs replicate in large reference?
    2. Cohort eQTL vs OneK1K — same, against single-cell reference
    3. Cohort eQTL vs AlphaGenome — does AlphaGenome predict cohort-observed effects?
    """
    log("\n" + "=" * 60)
    log("SANITY CHECK: Cohort-Level eQTLs vs Reference + AlphaGenome")
    log("=" * 60)

    sanity_metrics = {}

    # Create position key for cohort eQTLs (using hg38 coordinates)
    cohort_eqtls = cohort_eqtls.copy()
    cohort_eqtls['pos_key'] = cohort_eqtls['chrom'].astype(str) + '_' + cohort_eqtls['pos_hg38'].astype(int).astype(str)

    # --- 1. Cohort eQTLs vs DICE ---
    if dice_df is not None:
        log("\n1. Cohort eQTLs vs DICE...")
        dice_copy = dice_df.copy()
        dice_copy['pos_key'] = dice_copy['chrom'].astype(str) + '_' + dice_copy['pos'].astype(str)
        dice_best = dice_copy.loc[dice_copy.groupby('pos_key')['pval'].idxmin()]

        merged = cohort_eqtls.merge(
            dice_best[['pos_key', 'beta', 'dice_celltype']].rename(columns={'beta': 'dice_beta'}),
            on='pos_key', how='inner'
        )

        if len(merged) > 0:
            valid = merged[merged['beta'].notna() & merged['dice_beta'].notna()
                          & np.isfinite(merged['beta']) & np.isfinite(merged['dice_beta'])]
            if len(valid) > 0:
                concordance = (np.sign(valid['beta']) == np.sign(valid['dice_beta'])).mean()
                corr, pval = stats.spearmanr(valid['beta'].values, valid['dice_beta'].values, nan_policy='omit')

                sanity_metrics['cohort_vs_dice'] = {
                    'n_matched': len(valid),
                    'direction_concordance': float(concordance),
                    'spearman_r': float(corr) if not np.isnan(corr) else 0.0,
                    'spearman_pval': float(pval) if not np.isnan(pval) else 1.0,
                }
                log(f"  Matched: {len(valid):,}, concordance: {100*concordance:.1f}%, r={corr:.3f}")

                # Per-cohort breakdown
                per_cohort = {}
                for cohort in valid['cohort'].unique():
                    c = valid[valid['cohort'] == cohort]
                    if len(c) >= 10:
                        c_conc = (np.sign(c['beta']) == np.sign(c['dice_beta'])).mean()
                        per_cohort[cohort] = {'n': len(c), 'concordance': float(c_conc)}
                sanity_metrics['cohort_vs_dice']['per_cohort'] = per_cohort
            else:
                log("  No valid matched pairs after filtering")
        else:
            log("  No position matches found")

    # --- 2. Cohort eQTLs vs OneK1K ---
    if onek1k_df is not None:
        log("\n2. Cohort eQTLs vs OneK1K...")
        onek1k_copy = onek1k_df.copy()
        onek1k_copy['pos_key'] = onek1k_copy['chrom'].astype(str) + '_' + onek1k_copy['pos'].astype(str)
        onek1k_best = onek1k_copy.loc[onek1k_copy.groupby('pos_key')['pval'].idxmin()]

        merged = cohort_eqtls.merge(
            onek1k_best[['pos_key', 'beta', 'onek1k_celltype']].rename(columns={'beta': 'onek1k_beta'}),
            on='pos_key', how='inner'
        )

        if len(merged) > 0:
            valid = merged[merged['beta'].notna() & merged['onek1k_beta'].notna()
                          & np.isfinite(merged['beta']) & np.isfinite(merged['onek1k_beta'])]
            if len(valid) > 0:
                concordance = (np.sign(valid['beta']) == np.sign(valid['onek1k_beta'])).mean()
                corr, pval = stats.spearmanr(valid['beta'].values, valid['onek1k_beta'].values, nan_policy='omit')

                sanity_metrics['cohort_vs_onek1k'] = {
                    'n_matched': len(valid),
                    'direction_concordance': float(concordance),
                    'spearman_r': float(corr) if not np.isnan(corr) else 0.0,
                    'spearman_pval': float(pval) if not np.isnan(pval) else 1.0,
                }
                log(f"  Matched: {len(valid):,}, concordance: {100*concordance:.1f}%, r={corr:.3f}")

                # Per-cohort breakdown
                per_cohort = {}
                for cohort in valid['cohort'].unique():
                    c = valid[valid['cohort'] == cohort]
                    if len(c) >= 10:
                        c_conc = (np.sign(c['beta']) == np.sign(c['onek1k_beta'])).mean()
                        per_cohort[cohort] = {'n': len(c), 'concordance': float(c_conc)}
                sanity_metrics['cohort_vs_onek1k']['per_cohort'] = per_cohort
            else:
                log("  No valid matched pairs after filtering")
        else:
            log("  No position matches found")

    # --- 3. Cohort eQTLs vs AlphaGenome predictions ---
    log("\n3. Cohort eQTLs vs AlphaGenome predictions...")
    variants_copy = variants_df.copy()
    variants_copy['pos_key'] = variants_copy['chrom'].astype(str) + '_' + variants_copy['pos'].astype(str)

    # Select AlphaGenome impact column
    ag_col = None
    for candidate in ['ag_Whole_Blood', 'ag_Lymphocytes_EBV', 'alphagenome_impact_mean']:
        if candidate in variants_copy.columns:
            ag_col = candidate
            break

    if ag_col:
        merged = cohort_eqtls.merge(
            variants_copy[['pos_key', ag_col]],
            on='pos_key', how='inner'
        )

        if len(merged) > 0:
            # Filter to variants with non-zero AlphaGenome predictions
            valid = merged[merged[ag_col].notna() & (merged[ag_col] != 0)
                          & merged['beta'].notna() & np.isfinite(merged['beta'])]
            if len(valid) > 0:
                concordance = (np.sign(valid['beta']) == np.sign(valid[ag_col])).mean()
                corr, pval = stats.spearmanr(valid['beta'].values, valid[ag_col].values, nan_policy='omit')

                sanity_metrics['cohort_vs_alphagenome'] = {
                    'n_matched': len(valid),
                    'ag_column_used': ag_col,
                    'direction_concordance': float(concordance),
                    'spearman_r': float(corr) if not np.isnan(corr) else 0.0,
                    'spearman_pval': float(pval) if not np.isnan(pval) else 1.0,
                }
                log(f"  Matched: {len(valid):,}, concordance: {100*concordance:.1f}%, r={corr:.3f}")
            else:
                log("  No valid matched pairs with non-zero AlphaGenome predictions")
        else:
            log("  No position matches found")
    else:
        log("  No AlphaGenome impact column found")

    # --- Summary ---
    n_cohort_total = cohort_eqtls['pos_key'].nunique()
    sanity_metrics['summary'] = {
        'total_cohort_eqtl_positions': n_cohort_total,
        'total_cohort_eqtl_records': len(cohort_eqtls),
        'n_cohorts': int(cohort_eqtls['cohort'].nunique()),
    }

    return sanity_metrics


def compute_concordance(variants_df: pd.DataFrame) -> Dict:
    """
    Compute concordance metrics between AlphaGenome predictions and eQTL effects.
    """
    metrics = {}

    # DICE concordance
    dice_matched = variants_df[variants_df['dice_matched']].copy()
    if len(dice_matched) > 0:
        # Direction concordance: does AlphaGenome impact direction match DICE beta sign?
        if 'ag_Whole_Blood' in dice_matched.columns:
            ag_col = 'ag_Whole_Blood'
        elif 'ag_Lymphocytes_EBV' in dice_matched.columns:
            ag_col = 'ag_Lymphocytes_EBV'
        else:
            ag_col = 'alphagenome_impact_mean'

        same_direction = (np.sign(dice_matched[ag_col]) == np.sign(dice_matched['dice_beta']))
        concordance = same_direction.mean()

        # Correlation
        corr, pval = stats.spearmanr(
            dice_matched[ag_col].values,
            dice_matched['dice_beta'].values,
            nan_policy='omit'
        )

        metrics['dice'] = {
            'n_matched': len(dice_matched),
            'direction_concordance': float(concordance),
            'spearman_r': float(corr) if not np.isnan(corr) else 0.0,
            'spearman_pval': float(pval) if not np.isnan(pval) else 1.0,
        }
        log(f"  DICE: {len(dice_matched)} matches, {100*concordance:.1f}% concordance, r={corr:.3f}")

    # OneK1K concordance
    onek1k_matched = variants_df[variants_df['onek1k_matched']].copy()
    if len(onek1k_matched) > 0:
        if 'ag_Whole_Blood' in onek1k_matched.columns:
            ag_col = 'ag_Whole_Blood'
        elif 'ag_Lymphocytes_EBV' in onek1k_matched.columns:
            ag_col = 'ag_Lymphocytes_EBV'
        else:
            ag_col = 'alphagenome_impact_mean'

        same_direction = (np.sign(onek1k_matched[ag_col]) == np.sign(onek1k_matched['onek1k_beta']))
        concordance = same_direction.mean()

        corr, pval = stats.spearmanr(
            onek1k_matched[ag_col].values,
            onek1k_matched['onek1k_beta'].values,
            nan_policy='omit'
        )

        metrics['onek1k'] = {
            'n_matched': len(onek1k_matched),
            'direction_concordance': float(concordance),
            'spearman_r': float(corr) if not np.isnan(corr) else 0.0,
            'spearman_pval': float(pval) if not np.isnan(pval) else 1.0,
        }
        log(f"  OneK1K: {len(onek1k_matched)} matches, {100*concordance:.1f}% concordance, r={corr:.3f}")

    # Overall
    any_matched = variants_df['dice_matched'] | variants_df['onek1k_matched']
    metrics['overall'] = {
        'total_variants': len(variants_df),
        'any_eqtl_matched': int(any_matched.sum()),
        'both_matched': int((variants_df['dice_matched'] & variants_df['onek1k_matched']).sum()),
        'dice_only': int((variants_df['dice_matched'] & ~variants_df['onek1k_matched']).sum()),
        'onek1k_only': int((~variants_df['dice_matched'] & variants_df['onek1k_matched']).sum()),
    }

    return metrics


def generate_report(variants_df: pd.DataFrame, metrics: Dict, sanity_metrics: Optional[Dict], output_path: Path):
    """Generate markdown validation report."""
    log(f"Generating report: {output_path}")

    n_total = len(variants_df)
    n_dice = metrics.get('dice', {}).get('n_matched', 0)
    n_onek1k = metrics.get('onek1k', {}).get('n_matched', 0)

    report = f"""# eQTL Validation Report for Immunotherapy Variants
## Prioritizing Immunotherapy-Associated Regulatory Germline Variants Through AlphaGenome

## Validation Strategy

AlphaGenome is trained on GTEx RNA-seq data — validating against GTEx is circular.
Only **independent** eQTL databases are used for validation:

| Source | Independent? | Used? | Cell types |
|--------|-------------|-------|------------|
| DICE | Yes | Yes | 12 sorted immune cell types |
| OneK1K | Yes | Yes | 14 PBMC sc-eQTL cell types |
| GTEx | **No** (training data) | **Excluded** | — |

## Summary

| Metric | Value |
|--------|-------|
| Total variants | {n_total:,} |
| DICE matched | {n_dice:,} ({100*n_dice/n_total:.1f}%) |
| OneK1K matched | {n_onek1k:,} ({100*n_onek1k/n_total:.1f}%) |
| Either matched | {metrics['overall']['any_eqtl_matched']:,} ({100*metrics['overall']['any_eqtl_matched']/n_total:.1f}%) |
| Both matched | {metrics['overall']['both_matched']:,} |

## AlphaGenome vs DICE Concordance

"""

    if 'dice' in metrics:
        dice = metrics['dice']
        report += f"""| Metric | Value |
|--------|-------|
| Matched variants | {dice['n_matched']:,} |
| Direction concordance | {100*dice['direction_concordance']:.1f}% |
| Spearman correlation | {dice['spearman_r']:.3f} (p={dice['spearman_pval']:.2e}) |

"""

    report += """## AlphaGenome vs OneK1K Concordance

"""

    if 'onek1k' in metrics:
        onek1k = metrics['onek1k']
        report += f"""| Metric | Value |
|--------|-------|
| Matched variants | {onek1k['n_matched']:,} |
| Direction concordance | {100*onek1k['direction_concordance']:.1f}% |
| Spearman correlation | {onek1k['spearman_r']:.3f} (p={onek1k['spearman_pval']:.2e}) |

"""

    # Cohort eQTL sanity check section
    if sanity_metrics:
        report += """## Sanity Check: Cohort-Level eQTLs vs Reference Databases + AlphaGenome

Per-cohort eQTL analyses (from WES immunotherapy cohorts) are compared against
independent reference databases and AlphaGenome predictions. This validates that:
1. Small-cohort eQTLs replicate in large reference databases
2. AlphaGenome predictions agree with observed cohort-level eQTL effects

"""

        if 'cohort_vs_dice' in sanity_metrics:
            m = sanity_metrics['cohort_vs_dice']
            report += f"""### Cohort eQTLs vs DICE

| Metric | Value |
|--------|-------|
| Matched pairs | {m['n_matched']:,} |
| Direction concordance | {100*m['direction_concordance']:.1f}% |
| Spearman correlation | {m['spearman_r']:.3f} (p={m['spearman_pval']:.2e}) |

"""
            if 'per_cohort' in m and m['per_cohort']:
                report += "**Per-cohort breakdown (n >= 10):**\n\n"
                report += "| Cohort | N matched | Concordance |\n|--------|-----------|-------------|\n"
                for cohort, vals in sorted(m['per_cohort'].items(), key=lambda x: -x[1]['concordance']):
                    report += f"| {cohort} | {vals['n']:,} | {100*vals['concordance']:.1f}% |\n"
                report += "\n"

        if 'cohort_vs_onek1k' in sanity_metrics:
            m = sanity_metrics['cohort_vs_onek1k']
            report += f"""### Cohort eQTLs vs OneK1K

| Metric | Value |
|--------|-------|
| Matched pairs | {m['n_matched']:,} |
| Direction concordance | {100*m['direction_concordance']:.1f}% |
| Spearman correlation | {m['spearman_r']:.3f} (p={m['spearman_pval']:.2e}) |

"""

        if 'cohort_vs_alphagenome' in sanity_metrics:
            m = sanity_metrics['cohort_vs_alphagenome']
            report += f"""### Cohort eQTLs vs AlphaGenome Predictions

| Metric | Value |
|--------|-------|
| Matched pairs | {m['n_matched']:,} |
| AlphaGenome column | {m['ag_column_used']} |
| Direction concordance | {100*m['direction_concordance']:.1f}% |
| Spearman correlation | {m['spearman_r']:.3f} (p={m['spearman_pval']:.2e}) |

"""

    # Top variants by cohort
    report += """## Top Variants by Cohort

Showing variants with highest AlphaGenome impact that are also validated in eQTL databases:

"""

    any_matched = variants_df['dice_matched'] | variants_df['onek1k_matched']
    validated = variants_df[any_matched].sort_values('alphagenome_impact_mean', ascending=False)

    if len(validated) > 0:
        for cohort in validated['cohort'].unique()[:5]:
            cohort_vars = validated[validated['cohort'] == cohort].head(5)
            report += f"\n### {cohort}\n\n"
            report += "| rsid | gene | pval | AG impact | DICE | OneK1K |\n"
            report += "|------|------|------|-----------|------|--------|\n"

            for _, row in cohort_vars.iterrows():
                dice_status = f"✓ {row.get('dice_gene', '')}" if row.get('dice_matched', False) else "-"
                onek1k_status = f"✓ {row.get('onek1k_gene', '')}" if row.get('onek1k_matched', False) else "-"
                report += f"| {row.get('rsid', 'NA')} | {row.get('gene', 'NA')} | {row.get('pval', 0):.2e} | {row.get('alphagenome_impact_mean', 0):.4f} | {dice_status} | {onek1k_status} |\n"

    report += f"""

## Interpretation

This analysis validates AlphaGenome regulatory predictions against independent eQTL databases.
GTEx is excluded because AlphaGenome is trained on GTEx RNA-seq data (circular validation).

**Key findings:**
- {n_dice:,} therapy-associated variants ({100*n_dice/n_total:.1f}%) are known DICE immune cell eQTLs
- {n_onek1k:,} therapy-associated variants ({100*n_onek1k/n_total:.1f}%) are known OneK1K sc-eQTLs

**Note on AlphaGenome use case:**
This is the CORRECT use of AlphaGenome — predicting regulatory effects for variants with
UNKNOWN mechanisms (therapy association ≠ regulatory evidence). Matching to independent eQTL
databases (DICE, OneK1K) provides validation that AlphaGenome-predicted regulatory variants
have measurable eQTL effects in immune cells.

---
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(output_path, 'w') as f:
        f.write(report)


def main():
    parser = argparse.ArgumentParser(description='Validate AlphaGenome predictions against eQTL databases')
    parser.add_argument('--input', type=str, default=None,
                       help='Input h5ad file (default: alphagenome_predictions.h5ad)')
    parser.add_argument('--skip-onek1k', action='store_true',
                       help='Skip OneK1K matching')
    parser.add_argument('--skip-cohort-sanity', action='store_true',
                       help='Skip cohort-level eQTL sanity check')
    args = parser.parse_args()

    log("=" * 60)
    log("EQTL VALIDATION: IMMUNOTHERAPY VARIANTS")
    log("GTEx EXCLUDED (circular with AlphaGenome training)")
    log("Validation sources: DICE + OneK1K (independent)")
    log("=" * 60)

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load AlphaGenome predictions
    input_path = Path(args.input) if args.input else ALPHAGENOME_DIR / 'alphagenome_predictions.h5ad'
    variants_df = load_alphagenome_predictions(input_path)

    # Load eQTL databases
    log("\nLoading eQTL databases...")

    dice_df = load_dice_eqtls()
    if dice_df is not None:
        variants_df = match_to_dice(variants_df, dice_df)
    else:
        variants_df['dice_matched'] = False

    onek1k_df = None
    if not args.skip_onek1k:
        onek1k_df = load_onek1k_eqtls()
        if onek1k_df is not None:
            variants_df = match_to_onek1k(variants_df, onek1k_df)
        else:
            variants_df['onek1k_matched'] = False
    else:
        variants_df['onek1k_matched'] = False

    # Compute concordance metrics
    log("\nComputing concordance metrics...")
    metrics = compute_concordance(variants_df)

    # Cohort-level eQTL sanity check
    sanity_metrics = None
    if not args.skip_cohort_sanity:
        cohort_eqtls = load_cohort_eqtls()
        if cohort_eqtls is not None:
            sanity_metrics = run_cohort_sanity_check(
                variants_df, cohort_eqtls, dice_df, onek1k_df
            )

    # Save results
    log("\nSaving results...")

    # Save matched variants
    output_csv = OUTPUT_DIR / 'eqtl_matched.csv'
    variants_df.to_csv(output_csv, index=False)
    log(f"  Saved: {output_csv}")

    # Save metrics
    all_metrics = {'validation': metrics}
    if sanity_metrics:
        all_metrics['cohort_sanity_check'] = sanity_metrics

    metrics_path = OUTPUT_DIR / 'validation_metrics.json'
    with open(metrics_path, 'w') as f:
        json.dump(all_metrics, f, indent=2)
    log(f"  Saved: {metrics_path}")

    # Generate report
    report_path = OUTPUT_DIR / 'validation_report.md'
    generate_report(variants_df, metrics, sanity_metrics, report_path)
    log(f"  Saved: {report_path}")

    # Summary
    log("\n" + "=" * 60)
    log("COMPLETE")
    log("=" * 60)
    log(f"  Total variants: {len(variants_df):,}")
    log(f"  DICE matched: {metrics.get('dice', {}).get('n_matched', 0):,}")
    log(f"  OneK1K matched: {metrics.get('onek1k', {}).get('n_matched', 0):,}")
    if sanity_metrics:
        log(f"  Cohort sanity check: {sanity_metrics.get('summary', {}).get('n_cohorts', 0)} cohorts")

    log("\nNext steps:")
    log("  1. Run 04_prioritize_variants.py to compute final priority scores")


if __name__ == '__main__':
    main()
