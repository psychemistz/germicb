#!/usr/bin/env python3
"""
Validate AlphaGenome Predictions Against eQTL Databases
========================================================
Match immunotherapy-associated variants against multiple eQTL sources:

1. DICE - Immune cell type-specific eQTLs (12 cell types, hg38)
2. GTEx - Whole blood bulk eQTLs (v10, hg38)
3. CIMA - Single-cell immune eQTLs (69 cell types, hg38)

This validates AlphaGenome predictions for variants with unknown regulatory effects
against measured eQTL databases.

Input:
- results/alphagenome/alphagenome_predictions.h5ad
- DICE hg38 data (from cytokine atlas liftover)
- GTEx v10 whole blood

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

# eQTL database paths (from cytokine atlas project)
CYTOKINE_ATLAS_DIR = Path('/vf/users/parks34/projects/2secactpy/results/alphagenome')

# DICE hg38 (lifted from hg19)
DICE_HG38_DIR = CYTOKINE_ATLAS_DIR / 'dice_data' / 'hg38'

# GTEx v10 whole blood
GTEX_DIR = Path('/data/parks34/projects/2secactpy/results/alphagenome/gtex_data')
GTEX_PARQUET = GTEX_DIR / 'GTEx_Analysis_v10_eQTL_updated' / 'Whole_Blood.v10.eQTLs.signif_pairs.parquet'

# CIMA eQTLs (optional)
CIMA_EQTL_DIR = Path('/data/Jiang_Lab/Data/Seongyong/CIMA/xQTL')


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


def load_gtex_eqtls() -> Optional[pd.DataFrame]:
    """
    Load GTEx v10 whole blood eQTLs from parquet file.

    Returns DataFrame with columns: chrom, pos, ref, alt, gene_id, slope, pval_nominal
    """
    if not GTEX_PARQUET.exists():
        log(f"  GTEx parquet file not found: {GTEX_PARQUET}")
        return None

    log(f"Loading GTEx eQTLs: {GTEX_PARQUET}")

    gtex_df = pd.read_parquet(GTEX_PARQUET)
    log(f"  Loaded {len(gtex_df):,} GTEx eQTLs")

    # Parse variant_id: chr1_123456_A_G_b38
    parsed = gtex_df['variant_id'].str.extract(r'^(chr\d+)_(\d+)_([ACGT]+)_([ACGT]+)_b38$')
    gtex_df['chrom'] = parsed[0]
    gtex_df['pos'] = pd.to_numeric(parsed[1], errors='coerce')
    gtex_df['ref'] = parsed[2]
    gtex_df['alt'] = parsed[3]

    # Filter to valid coordinates
    gtex_df = gtex_df[gtex_df['chrom'].notna()].copy()
    log(f"  Valid coordinates: {len(gtex_df):,}")

    return gtex_df


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


def match_to_gtex(variants_df: pd.DataFrame, gtex_df: pd.DataFrame) -> pd.DataFrame:
    """
    Match therapy variants to GTEx eQTLs by position.

    Returns variants_df with added GTEx columns.
    """
    log("Matching to GTEx eQTLs...")

    # Create position keys - ensure integer positions for consistent matching
    variants_df = variants_df.copy()
    variants_df['pos_key'] = variants_df['chrom'].astype(str) + '_' + variants_df['pos'].astype(int).astype(str)

    gtex_df = gtex_df.copy()
    gtex_df['pos_key'] = gtex_df['chrom'].astype(str) + '_' + gtex_df['pos'].astype(int).astype(str)

    # Find matches
    common_positions = set(variants_df['pos_key']) & set(gtex_df['pos_key'])
    log(f"  Matching positions: {len(common_positions)}")

    if len(common_positions) == 0:
        variants_df['gtex_matched'] = False
        return variants_df

    # Get best GTEx hit per position (lowest p-value)
    gtex_best = gtex_df.loc[gtex_df.groupby('pos_key')['pval_nominal'].idxmin()]

    # Merge
    gtex_subset = gtex_best[['pos_key', 'gene_id', 'slope', 'pval_nominal']].copy()
    gtex_subset = gtex_subset.rename(columns={
        'gene_id': 'gtex_gene',
        'slope': 'gtex_slope',
        'pval_nominal': 'gtex_pval',
    })

    variants_df = variants_df.merge(gtex_subset, on='pos_key', how='left')
    variants_df['gtex_matched'] = variants_df['gtex_gene'].notna()

    n_matched = variants_df['gtex_matched'].sum()
    log(f"  Matched: {n_matched:,} / {len(variants_df):,} ({100*n_matched/len(variants_df):.1f}%)")

    return variants_df


def compute_concordance(variants_df: pd.DataFrame) -> Dict:
    """
    Compute concordance metrics between AlphaGenome predictions and eQTL effects.
    """
    metrics = {}

    # DICE concordance
    dice_matched = variants_df[variants_df['dice_matched']].copy()
    if len(dice_matched) > 0:
        # Direction concordance: does AlphaGenome impact direction match DICE beta sign?
        # For RNA-seq tracks, positive diff = increased expression
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

    # GTEx concordance
    gtex_matched = variants_df[variants_df['gtex_matched']].copy()
    if len(gtex_matched) > 0:
        if 'ag_Whole_Blood' in gtex_matched.columns:
            ag_col = 'ag_Whole_Blood'
        elif 'ag_Lymphocytes_EBV' in gtex_matched.columns:
            ag_col = 'ag_Lymphocytes_EBV'
        else:
            ag_col = 'alphagenome_impact_mean'

        same_direction = (np.sign(gtex_matched[ag_col]) == np.sign(gtex_matched['gtex_slope']))
        concordance = same_direction.mean()

        corr, pval = stats.spearmanr(
            gtex_matched[ag_col].values,
            gtex_matched['gtex_slope'].values,
            nan_policy='omit'
        )

        metrics['gtex'] = {
            'n_matched': len(gtex_matched),
            'direction_concordance': float(concordance),
            'spearman_r': float(corr) if not np.isnan(corr) else 0.0,
            'spearman_pval': float(pval) if not np.isnan(pval) else 1.0,
        }
        log(f"  GTEx: {len(gtex_matched)} matches, {100*concordance:.1f}% concordance, r={corr:.3f}")

    # Overall
    any_matched = variants_df['dice_matched'] | variants_df['gtex_matched']
    metrics['overall'] = {
        'total_variants': len(variants_df),
        'any_eqtl_matched': int(any_matched.sum()),
        'both_matched': int((variants_df['dice_matched'] & variants_df['gtex_matched']).sum()),
        'dice_only': int((variants_df['dice_matched'] & ~variants_df['gtex_matched']).sum()),
        'gtex_only': int((~variants_df['dice_matched'] & variants_df['gtex_matched']).sum()),
    }

    return metrics


def generate_report(variants_df: pd.DataFrame, metrics: Dict, output_path: Path):
    """Generate markdown validation report."""
    log(f"Generating report: {output_path}")

    n_total = len(variants_df)
    n_dice = metrics.get('dice', {}).get('n_matched', 0)
    n_gtex = metrics.get('gtex', {}).get('n_matched', 0)

    report = f"""# eQTL Validation Report for Immunotherapy Variants

## Summary

| Metric | Value |
|--------|-------|
| Total variants | {n_total:,} |
| DICE matched | {n_dice:,} ({100*n_dice/n_total:.1f}%) |
| GTEx matched | {n_gtex:,} ({100*n_gtex/n_total:.1f}%) |
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

    report += """## AlphaGenome vs GTEx Concordance

"""

    if 'gtex' in metrics:
        gtex = metrics['gtex']
        report += f"""| Metric | Value |
|--------|-------|
| Matched variants | {gtex['n_matched']:,} |
| Direction concordance | {100*gtex['direction_concordance']:.1f}% |
| Spearman correlation | {gtex['spearman_r']:.3f} (p={gtex['spearman_pval']:.2e}) |

"""

    # Top variants by cohort
    report += """## Top Variants by Cohort

Showing variants with highest AlphaGenome impact that are also validated in eQTL databases:

"""

    any_matched = variants_df['dice_matched'] | variants_df['gtex_matched']
    validated = variants_df[any_matched].sort_values('alphagenome_impact_mean', ascending=False)

    if len(validated) > 0:
        for cohort in validated['cohort'].unique()[:5]:
            cohort_vars = validated[validated['cohort'] == cohort].head(5)
            report += f"\n### {cohort}\n\n"
            report += "| rsid | gene | pval | AG impact | DICE | GTEx |\n"
            report += "|------|------|------|-----------|------|------|\n"

            for _, row in cohort_vars.iterrows():
                dice_status = f"✓ {row.get('dice_gene', '')}" if row.get('dice_matched', False) else "-"
                gtex_status = f"✓ {row.get('gtex_gene', '')[:15]}" if row.get('gtex_matched', False) else "-"
                report += f"| {row.get('rsid', 'NA')} | {row.get('gene', 'NA')} | {row.get('pval', 0):.2e} | {row.get('alphagenome_impact_mean', 0):.4f} | {dice_status} | {gtex_status} |\n"

    report += f"""

## Interpretation

This analysis validates AlphaGenome regulatory predictions against measured eQTL databases.

**Key findings:**
- {n_dice:,} therapy-associated variants ({100*n_dice/n_total:.1f}%) are known DICE immune cell eQTLs
- {n_gtex:,} therapy-associated variants ({100*n_gtex/n_total:.1f}%) are known GTEx whole blood eQTLs

**Note on AlphaGenome use case:**
This is the CORRECT use of AlphaGenome - predicting regulatory effects for variants with
UNKNOWN mechanisms (therapy association ≠ regulatory evidence). Matching to eQTL databases
provides independent validation that AlphaGenome-predicted regulatory variants have measurable
eQTL effects in immune cells.

---
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(output_path, 'w') as f:
        f.write(report)


def main():
    parser = argparse.ArgumentParser(description='Validate AlphaGenome predictions against eQTL databases')
    parser.add_argument('--input', type=str, default=None,
                       help='Input h5ad file (default: alphagenome_predictions.h5ad)')
    parser.add_argument('--skip-gtex', action='store_true',
                       help='Skip GTEx matching (faster)')
    args = parser.parse_args()

    log("=" * 60)
    log("EQTL VALIDATION: IMMUNOTHERAPY VARIANTS")
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

    if not args.skip_gtex:
        gtex_df = load_gtex_eqtls()
        if gtex_df is not None:
            variants_df = match_to_gtex(variants_df, gtex_df)
        else:
            variants_df['gtex_matched'] = False
    else:
        variants_df['gtex_matched'] = False

    # Compute concordance metrics
    log("\nComputing concordance metrics...")
    metrics = compute_concordance(variants_df)

    # Save results
    log("\nSaving results...")

    # Save matched variants
    output_csv = OUTPUT_DIR / 'eqtl_matched.csv'
    variants_df.to_csv(output_csv, index=False)
    log(f"  Saved: {output_csv}")

    # Save metrics
    metrics_path = OUTPUT_DIR / 'validation_metrics.json'
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)
    log(f"  Saved: {metrics_path}")

    # Generate report
    report_path = OUTPUT_DIR / 'validation_report.md'
    generate_report(variants_df, metrics, report_path)
    log(f"  Saved: {report_path}")

    # Summary
    log("\n" + "=" * 60)
    log("COMPLETE")
    log("=" * 60)
    log(f"  Total variants: {len(variants_df):,}")
    log(f"  DICE matched: {metrics.get('dice', {}).get('n_matched', 0):,}")
    log(f"  GTEx matched: {metrics.get('gtex', {}).get('n_matched', 0):,}")

    log("\nNext steps:")
    log("  1. Run 04_prioritize_variants.py to compute final priority scores")


if __name__ == '__main__':
    main()
