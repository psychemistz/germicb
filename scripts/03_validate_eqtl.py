#!/usr/bin/env python3
"""
Validate AlphaGenome Predictions Against eQTL Databases
========================================================
Prioritizing Immunotherapy-Associated Regulatory Germline Variants Through AlphaGenome

Match immunotherapy-associated variants against independent eQTL sources:

1. DICE - Immune cell type-specific eQTLs (12 cell types, hg38)
2. OneK1K - Single-cell PBMC eQTLs (14 cell types, hg38; Yazar et al., Science 2022)
3. Cohort eQTLs - Per-cohort local eQTL analysis (sanity check, hg19->hg38 liftover)

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

import sys
import json
import time
import argparse
from pathlib import Path
from typing import Dict, Optional

# Ensure project root is on sys.path for lib imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy import stats
import anndata as ad

from lib.config import (
    ALPHAGENOME_DIR, EQTL_VALIDATION_DIR as OUTPUT_DIR,
)
from lib.log import log
from lib.variants import make_pos_key
from lib.eqtl import load_dice, load_onek1k, load_cohort_eqtls
from lib.matching import match_by_position, compute_concordance


def load_alphagenome_predictions(input_path: Path) -> pd.DataFrame:
    """
    Load AlphaGenome predictions from h5ad file.

    Returns DataFrame with variant info and track predictions.
    """
    log(f"Loading AlphaGenome predictions: {input_path}")

    adata = ad.read_h5ad(input_path)
    log(f"  Loaded {adata.n_obs:,} variants, {adata.n_vars} tracks")

    variants_df = adata.obs.copy()
    X = adata.X

    variants_df['alphagenome_impact_mean'] = np.abs(X).mean(axis=1)
    variants_df['alphagenome_impact_max'] = np.abs(X).max(axis=1)

    track_names = adata.var_names.tolist()
    for i, track in enumerate(track_names):
        clean_name = track.replace('RNA_SEQ_GTEx_', '').replace(' ', '_')
        variants_df[f'ag_{clean_name}'] = X[:, i]

    log(f"  Impact score range: {variants_df['alphagenome_impact_mean'].min():.4f} to {variants_df['alphagenome_impact_mean'].max():.4f}")

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

    cohort_eqtls = cohort_eqtls.copy()
    cohort_eqtls['pos_key'] = cohort_eqtls['chrom'].astype(str) + '_' + cohort_eqtls['pos_hg38'].astype(int).astype(str)

    # --- 1. Cohort eQTLs vs DICE ---
    if dice_df is not None:
        log("\n1. Cohort eQTLs vs DICE...")
        dice_copy = dice_df.copy()
        dice_copy['pos_key'] = make_pos_key(dice_copy)
        dice_best = dice_copy.loc[dice_copy.groupby('pos_key')['pval'].idxmin()]

        merged = cohort_eqtls.merge(
            dice_best[['pos_key', 'beta', 'dice_celltype']].rename(columns={'beta': 'dice_beta'}),
            on='pos_key', how='inner'
        )

        if len(merged) > 0:
            conc = compute_concordance(merged, 'beta', 'dice_beta')
            if conc['n_valid'] > 0:
                sanity_metrics['cohort_vs_dice'] = {
                    'n_matched': conc['n_valid'],
                    'direction_concordance': conc['direction_concordance'],
                    'spearman_r': conc['spearman_r'],
                    'spearman_pval': conc['spearman_pval'],
                }
                log(f"  Matched: {conc['n_valid']:,}, concordance: {100*conc['direction_concordance']:.1f}%, r={conc['spearman_r']:.3f}")

                # Per-cohort breakdown
                valid = merged[merged['beta'].notna() & merged['dice_beta'].notna()
                              & np.isfinite(merged['beta']) & np.isfinite(merged['dice_beta'])]
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
        onek1k_copy['pos_key'] = make_pos_key(onek1k_copy)
        onek1k_best = onek1k_copy.loc[onek1k_copy.groupby('pos_key')['pval'].idxmin()]

        merged = cohort_eqtls.merge(
            onek1k_best[['pos_key', 'beta', 'onek1k_celltype']].rename(columns={'beta': 'onek1k_beta'}),
            on='pos_key', how='inner'
        )

        if len(merged) > 0:
            conc = compute_concordance(merged, 'beta', 'onek1k_beta')
            if conc['n_valid'] > 0:
                sanity_metrics['cohort_vs_onek1k'] = {
                    'n_matched': conc['n_valid'],
                    'direction_concordance': conc['direction_concordance'],
                    'spearman_r': conc['spearman_r'],
                    'spearman_pval': conc['spearman_pval'],
                }
                log(f"  Matched: {conc['n_valid']:,}, concordance: {100*conc['direction_concordance']:.1f}%, r={conc['spearman_r']:.3f}")

                valid = merged[merged['beta'].notna() & merged['onek1k_beta'].notna()
                              & np.isfinite(merged['beta']) & np.isfinite(merged['onek1k_beta'])]
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
    variants_copy['pos_key'] = make_pos_key(variants_copy)

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
            conc = compute_concordance(merged, 'beta', ag_col)
            if conc['n_valid'] > 0:
                sanity_metrics['cohort_vs_alphagenome'] = {
                    'n_matched': conc['n_valid'],
                    'ag_column_used': ag_col,
                    'direction_concordance': conc['direction_concordance'],
                    'spearman_r': conc['spearman_r'],
                    'spearman_pval': conc['spearman_pval'],
                }
                log(f"  Matched: {conc['n_valid']:,}, concordance: {100*conc['direction_concordance']:.1f}%, r={conc['spearman_r']:.3f}")
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


def compute_validation_metrics(variants_df: pd.DataFrame) -> Dict:
    """
    Compute concordance metrics between AlphaGenome predictions and eQTL effects.
    """
    metrics = {}

    # DICE concordance
    dice_matched = variants_df[variants_df['dice_matched']].copy()
    if len(dice_matched) > 0:
        for ag_col in ['ag_Whole_Blood', 'ag_Lymphocytes_EBV', 'alphagenome_impact_mean']:
            if ag_col in dice_matched.columns:
                break

        conc = compute_concordance(dice_matched, ag_col, 'dice_beta')
        metrics['dice'] = {
            'n_matched': len(dice_matched),
            'direction_concordance': conc['direction_concordance'],
            'spearman_r': conc['spearman_r'],
            'spearman_pval': conc['spearman_pval'],
        }
        log(f"  DICE: {len(dice_matched)} matches, {100*conc['direction_concordance']:.1f}% concordance, r={conc['spearman_r']:.3f}")

    # OneK1K concordance
    onek1k_matched = variants_df[variants_df['onek1k_matched']].copy()
    if len(onek1k_matched) > 0:
        for ag_col in ['ag_Whole_Blood', 'ag_Lymphocytes_EBV', 'alphagenome_impact_mean']:
            if ag_col in onek1k_matched.columns:
                break

        conc = compute_concordance(onek1k_matched, ag_col, 'onek1k_beta')
        metrics['onek1k'] = {
            'n_matched': len(onek1k_matched),
            'direction_concordance': conc['direction_concordance'],
            'spearman_r': conc['spearman_r'],
            'spearman_pval': conc['spearman_pval'],
        }
        log(f"  OneK1K: {len(onek1k_matched)} matches, {100*conc['direction_concordance']:.1f}% concordance, r={conc['spearman_r']:.3f}")

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
                dice_status = f"+ {row.get('dice_gene', '')}" if row.get('dice_matched', False) else "-"
                onek1k_status = f"+ {row.get('onek1k_gene', '')}" if row.get('onek1k_matched', False) else "-"
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
UNKNOWN mechanisms (therapy association != regulatory evidence). Matching to independent eQTL
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

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load AlphaGenome predictions
    input_path = Path(args.input) if args.input else ALPHAGENOME_DIR / 'alphagenome_predictions.h5ad'
    variants_df = load_alphagenome_predictions(input_path)

    # Load eQTL databases and match
    log("\nLoading eQTL databases...")

    dice_df = load_dice()
    if dice_df is not None:
        variants_df = match_by_position(variants_df, dice_df, 'dice')
    else:
        variants_df['dice_matched'] = False

    onek1k_df = None
    if not args.skip_onek1k:
        onek1k_df = load_onek1k()
        if onek1k_df is not None:
            variants_df = match_by_position(variants_df, onek1k_df, 'onek1k')
        else:
            variants_df['onek1k_matched'] = False
    else:
        variants_df['onek1k_matched'] = False

    # Compute concordance metrics
    log("\nComputing concordance metrics...")
    metrics = compute_validation_metrics(variants_df)

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

    output_csv = OUTPUT_DIR / 'eqtl_matched.csv'
    variants_df.to_csv(output_csv, index=False)
    log(f"  Saved: {output_csv}")

    all_metrics = {'validation': metrics}
    if sanity_metrics:
        all_metrics['cohort_sanity_check'] = sanity_metrics

    metrics_path = OUTPUT_DIR / 'validation_metrics.json'
    with open(metrics_path, 'w') as f:
        json.dump(all_metrics, f, indent=2)
    log(f"  Saved: {metrics_path}")

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
