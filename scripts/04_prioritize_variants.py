#!/usr/bin/env python3
"""
Prioritize Immunotherapy-Associated Variants
=============================================
Compute composite priority scores using:

1. Therapy association strength (-log10 p-value)
2. AlphaGenome predicted regulatory impact
3. eQTL validation (DICE, OneK1K)

Tier Assignment:
- Tier 1: Therapy-associated + AlphaGenome-predicted + eQTL-validated
- Tier 2: Therapy-associated + AlphaGenome-predicted (no eQTL match)
- Tier 3: Therapy-associated + eQTL-validated (low AlphaGenome prediction)
- Tier 4: Therapy-associated only

Input:
- results/eqtl_validation/eqtl_matched.csv

Output:
- results/prioritized/prioritized_variants.csv
- results/prioritized/tier_summary.csv
- results/prioritized/prioritization_report.md
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

# ==============================================================================
# Configuration
# ==============================================================================

PROJECT_DIR = Path('/data/parks34/projects/4germicb')
INPUT_DIR = PROJECT_DIR / 'results' / 'eqtl_validation'
OUTPUT_DIR = PROJECT_DIR / 'results' / 'prioritized'

# Priority score thresholds
ALPHAGENOME_HIGH_THRESHOLD = 0.5      # Impact score for "high" prediction
ALPHAGENOME_MEDIUM_THRESHOLD = 0.2    # Impact score for "medium" prediction
PVAL_STRINGENT = 0.001                # Stringent p-value threshold


def log(msg: str):
    """Print timestamped log message."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def compute_priority_score(row: pd.Series) -> float:
    """
    Compute composite priority score.

    Components:
    - therapy_score: -log10(pval), capped at 10
    - alphagenome_score: normalized impact score (0-3)
    - eqtl_score: bonus for eQTL validation (0-2)

    Returns:
        Priority score (higher = higher priority)
    """
    score = 0.0

    # Therapy association score (-log10 p-value, capped)
    pval = row.get('pval', 1)
    if pval > 0:
        therapy_score = min(-np.log10(pval), 10)
    else:
        therapy_score = 10
    score += therapy_score

    # AlphaGenome impact score
    ag_impact = row.get('alphagenome_impact_mean', 0)
    if ag_impact >= ALPHAGENOME_HIGH_THRESHOLD:
        score += 3
    elif ag_impact >= ALPHAGENOME_MEDIUM_THRESHOLD:
        score += 2
    elif ag_impact > 0:
        score += 1

    # eQTL validation bonus
    if row.get('dice_matched', False):
        score += 2
    if row.get('onek1k_matched', False):
        score += 1

    return score


def assign_tier(row: pd.Series) -> int:
    """
    Assign evidence tier based on available evidence.

    Tiers:
    - 1: Therapy + AlphaGenome (high) + eQTL validated
    - 2: Therapy + AlphaGenome (high/medium) only
    - 3: Therapy + eQTL validated (low/no AlphaGenome)
    - 4: Therapy only
    """
    ag_impact = row.get('alphagenome_impact_mean', 0)
    has_eqtl = row.get('dice_matched', False) or row.get('onek1k_matched', False)
    has_high_ag = ag_impact >= ALPHAGENOME_MEDIUM_THRESHOLD

    if has_high_ag and has_eqtl:
        return 1
    elif has_high_ag:
        return 2
    elif has_eqtl:
        return 3
    else:
        return 4


def generate_report(df: pd.DataFrame, output_path: Path):
    """Generate markdown prioritization report."""
    log(f"Generating report: {output_path}")

    n_total = len(df)

    # Tier counts
    tier_counts = df['tier'].value_counts().sort_index()

    # Cohort breakdown
    cohort_summary = df.groupby('cohort').agg({
        'tier': lambda x: (x == 1).sum(),
        'priority_score': 'mean',
        'rsid': 'count'
    }).rename(columns={'tier': 'tier1_count', 'rsid': 'total'})

    report = f"""# Immunotherapy Variant Prioritization Report

## Summary

| Metric | Value |
|--------|-------|
| Total variants | {n_total:,} |
| Tier 1 (AG + eQTL) | {tier_counts.get(1, 0):,} ({100*tier_counts.get(1, 0)/n_total:.1f}%) |
| Tier 2 (AG only) | {tier_counts.get(2, 0):,} ({100*tier_counts.get(2, 0)/n_total:.1f}%) |
| Tier 3 (eQTL only) | {tier_counts.get(3, 0):,} ({100*tier_counts.get(3, 0)/n_total:.1f}%) |
| Tier 4 (therapy only) | {tier_counts.get(4, 0):,} ({100*tier_counts.get(4, 0)/n_total:.1f}%) |

## Tier Definitions

| Tier | Criteria | Interpretation |
|------|----------|----------------|
| 1 | Therapy + AlphaGenome + eQTL | Strong regulatory evidence, validated mechanism |
| 2 | Therapy + AlphaGenome | Predicted regulatory, not yet validated in eQTL databases |
| 3 | Therapy + eQTL | Known eQTL but AlphaGenome didn't predict high impact |
| 4 | Therapy only | Association only, no regulatory evidence |

## Cohort Summary

| Cohort | Total | Tier 1 | Mean Priority |
|--------|-------|--------|---------------|
"""

    for cohort, row in cohort_summary.iterrows():
        report += f"| {cohort} | {int(row['total']):,} | {int(row['tier1_count'])} | {row['priority_score']:.2f} |\n"

    # Top Tier 1 variants
    tier1 = df[df['tier'] == 1].nlargest(20, 'priority_score')

    report += f"""

## Top 20 Tier 1 Variants

These variants have the strongest evidence: therapy association + AlphaGenome prediction + eQTL validation.

| rsid | gene | cohort | pval | AG impact | DICE | OneK1K | Priority |
|------|------|--------|------|-----------|------|--------|----------|
"""

    for _, row in tier1.iterrows():
        dice_gene = row.get('dice_gene', '') if pd.notna(row.get('dice_gene', '')) else '-'
        onek1k_gene = str(row.get('onek1k_gene', ''))[:12] if pd.notna(row.get('onek1k_gene', '')) else '-'
        report += f"| {row['rsid']} | {row['gene']} | {row['cohort'][:20]} | {row['pval']:.2e} | {row['alphagenome_impact_mean']:.3f} | {dice_gene} | {onek1k_gene} | {row['priority_score']:.1f} |\n"

    # Top variants by cohort
    report += """

## Top Variants by Cohort

"""

    for cohort in df['cohort'].unique():
        cohort_df = df[df['cohort'] == cohort]
        top5 = cohort_df.nlargest(5, 'priority_score')

        report += f"\n### {cohort}\n\n"
        report += "| rsid | gene | pval | AG impact | Tier | Priority |\n"
        report += "|------|------|------|-----------|------|----------|\n"

        for _, row in top5.iterrows():
            report += f"| {row['rsid']} | {row['gene']} | {row['pval']:.2e} | {row['alphagenome_impact_mean']:.3f} | {row['tier']} | {row['priority_score']:.1f} |\n"

    # Meta-analysis candidates (variants in multiple cohorts)
    report += """

## Cross-Cohort Variants

Variants appearing in multiple cohorts with Tier 1 or 2 evidence:

"""

    # Find variants by position appearing in multiple cohorts
    df['pos_id'] = df['chrom'] + ':' + df['pos'].astype(str)
    multi_cohort = df.groupby('pos_id').filter(lambda x: len(x) > 1)

    if len(multi_cohort) > 0:
        # Get top by mean priority across cohorts
        mc_summary = multi_cohort.groupby('pos_id').agg({
            'rsid': 'first',
            'gene': 'first',
            'cohort': lambda x: ', '.join(x.unique()[:3]),
            'tier': 'min',
            'priority_score': 'mean',
        }).nlargest(15, 'priority_score')

        report += "| rsid | gene | cohorts | best tier | mean priority |\n"
        report += "|------|------|---------|-----------|---------------|\n"

        for pos_id, row in mc_summary.iterrows():
            report += f"| {row['rsid']} | {row['gene']} | {row['cohort']} | {row['tier']} | {row['priority_score']:.1f} |\n"
    else:
        report += "_No variants found in multiple cohorts_\n"

    report += f"""

## Next Steps

1. **Tier 1 variants**: Prioritize for functional validation
   - Gene expression studies in immune cells
   - Reporter assays for regulatory activity
   - CRISPR editing in relevant cell types

2. **Tier 2 variants**: Query real AlphaGenome API
   - Replace mock predictions with actual API calls
   - Re-evaluate with real regulatory predictions

3. **Cross-cohort variants**: Meta-analysis
   - Fixed-effects meta-analysis for consistent effect direction
   - Test replication across cancer types

---
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}
"""

    with open(output_path, 'w') as f:
        f.write(report)


def main():
    parser = argparse.ArgumentParser(description='Prioritize immunotherapy variants')
    parser.add_argument('--input', type=str, default=None,
                       help='Input CSV file (default: eqtl_matched.csv)')
    args = parser.parse_args()

    log("=" * 60)
    log("VARIANT PRIORITIZATION: IMMUNOTHERAPY COHORTS")
    log("=" * 60)

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load validated variants
    input_path = Path(args.input) if args.input else INPUT_DIR / 'eqtl_matched.csv'
    log(f"\nLoading variants: {input_path}")

    df = pd.read_csv(input_path)
    log(f"  Loaded {len(df):,} variants")

    # Compute priority scores
    log("\nComputing priority scores...")
    df['priority_score'] = df.apply(compute_priority_score, axis=1)
    df['tier'] = df.apply(assign_tier, axis=1)

    log(f"  Score range: {df['priority_score'].min():.1f} to {df['priority_score'].max():.1f}")

    # Tier breakdown
    log("\nTier breakdown:")
    for tier in [1, 2, 3, 4]:
        n = (df['tier'] == tier).sum()
        log(f"  Tier {tier}: {n:,} ({100*n/len(df):.1f}%)")

    # Sort by priority
    df = df.sort_values('priority_score', ascending=False)

    # Save prioritized variants
    output_csv = OUTPUT_DIR / 'prioritized_variants.csv'
    df.to_csv(output_csv, index=False)
    log(f"\nSaved: {output_csv}")

    # Save tier summary
    tier_summary = df.groupby(['cohort', 'tier']).size().unstack(fill_value=0)
    tier_summary_path = OUTPUT_DIR / 'tier_summary.csv'
    tier_summary.to_csv(tier_summary_path)
    log(f"Saved: {tier_summary_path}")

    # Generate report
    report_path = OUTPUT_DIR / 'prioritization_report.md'
    generate_report(df, report_path)
    log(f"Saved: {report_path}")

    # Summary
    log("\n" + "=" * 60)
    log("COMPLETE")
    log("=" * 60)

    # Show top 10 variants
    log("\nTop 10 prioritized variants:")
    log("-" * 80)
    top10 = df.head(10)
    for _, row in top10.iterrows():
        dice = "DICE" if row.get('dice_matched', False) else ""
        onek1k = "OneK1K" if row.get('onek1k_matched', False) else ""
        eqtl = f"[{dice}{'+' if dice and onek1k else ''}{onek1k}]" if dice or onek1k else ""
        log(f"  {row['rsid']:15} | {row['gene']:10} | {row['cohort'][:25]:25} | "
            f"p={row['pval']:.1e} | AG={row['alphagenome_impact_mean']:.2f} | "
            f"T{row['tier']} {eqtl}")

    log("-" * 80)
    log(f"\nTotal Tier 1 variants: {(df['tier'] == 1).sum()}")
    log(f"Report: {report_path}")


if __name__ == '__main__':
    main()
