#!/usr/bin/env python3
"""
Extract therapy-associated variants from WES immunotherapy cohorts.

This script aggregates variants associated with therapy response across
14 immunotherapy cohorts with consensus variant calling (GATK ∩ DeepVariant ∩ FreeBayes).

Input:
    - WES_ImmunoPredict/FINISHED/{cohort}/6logis_batch/output/fin_tbl.txt

Output:
    - data/therapy_variants_all.csv (all variants p < 0.05)
    - data/therapy_variants_stringent.csv (p < 0.001)
    - data/cohort_summary.csv (per-cohort statistics)
"""

import os
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd
import numpy as np

# ==============================================================================
# Configuration
# ==============================================================================

WES_DIR = Path('/data/parks34/projects/project_outdated/WES_ImmunoPredict/FINISHED')
OUTPUT_DIR = Path('/data/parks34/projects/4germicb/data')

# P-value thresholds
PVAL_SUGGESTIVE = 0.05
PVAL_STRINGENT = 0.001

# Directories to skip (not cohort directories)
SKIP_DIRS = {'summary_all', 'logis_batch', 'ridge', 'Transcriptome', 'DNAFormer', 'geneformer'}


def get_cohort_dirs() -> List[Path]:
    """Get list of cohort directories."""
    cohorts = []
    for d in WES_DIR.iterdir():
        if d.is_dir() and d.name not in SKIP_DIRS:
            # Check if it has the expected structure
            if (d / '6logis_batch').exists() or (d / '4merge_vcf').exists():
                cohorts.append(d)
    return sorted(cohorts)


def load_association_results(cohort_dir: Path) -> Optional[pd.DataFrame]:
    """Load association results from logis_batch output."""
    fin_tbl = cohort_dir / '6logis_batch' / 'output' / 'fin_tbl.txt'

    if not fin_tbl.exists():
        print(f"  [SKIP] No association results: {cohort_dir.name}")
        return None

    try:
        df = pd.read_csv(fin_tbl, sep='\t')
        df['cohort'] = cohort_dir.name

        # Standardize column names
        col_mapping = {
            'CHROM': 'chrom',
            'POS': 'pos',
            'REF': 'ref',
            'ALT': 'alt',
            'GENE': 'gene',
            'pvalue': 'pval',
            'coef': 'beta',
            'zscore': 'zscore'
        }
        df = df.rename(columns={k: v for k, v in col_mapping.items() if k in df.columns})

        return df

    except Exception as e:
        print(f"  [ERROR] Failed to load {fin_tbl}: {e}")
        return None


def load_variant_annotations(cohort_dir: Path) -> Optional[pd.DataFrame]:
    """Load variant annotations from merge_vcf output."""
    filt_var = cohort_dir / '4merge_vcf' / 'filt_var.txt'

    if not filt_var.exists():
        return None

    try:
        df = pd.read_csv(filt_var, sep='\t')
        return df
    except Exception as e:
        print(f"  [WARN] Failed to load annotations: {e}")
        return None


def main():
    print("=" * 60)
    print("Extract Therapy-Associated Variants from WES Cohorts")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cohorts = get_cohort_dirs()
    print(f"\nFound {len(cohorts)} cohorts")

    # Collect all results
    all_variants = []
    cohort_stats = []

    for cohort_dir in cohorts:
        print(f"\n[{cohort_dir.name}]")

        # Load association results
        df = load_association_results(cohort_dir)
        if df is None:
            continue

        n_total = len(df)
        n_suggestive = (df['pval'] < PVAL_SUGGESTIVE).sum() if 'pval' in df.columns else 0
        n_stringent = (df['pval'] < PVAL_STRINGENT).sum() if 'pval' in df.columns else 0

        print(f"  Total variants: {n_total:,}")
        print(f"  p < 0.05: {n_suggestive:,}")
        print(f"  p < 0.001: {n_stringent:,}")

        # Load annotations if available
        annot = load_variant_annotations(cohort_dir)
        if annot is not None:
            print(f"  Annotations available: {len(annot):,} variants")

        all_variants.append(df)
        cohort_stats.append({
            'cohort': cohort_dir.name,
            'n_total': n_total,
            'n_suggestive': n_suggestive,
            'n_stringent': n_stringent
        })

    if not all_variants:
        print("\nNo variants found!")
        return

    # Combine all variants
    combined = pd.concat(all_variants, ignore_index=True)
    print(f"\n{'=' * 60}")
    print(f"TOTAL: {len(combined):,} variants across {len(cohort_stats)} cohorts")

    # Filter by p-value
    if 'pval' in combined.columns:
        suggestive = combined[combined['pval'] < PVAL_SUGGESTIVE]
        stringent = combined[combined['pval'] < PVAL_STRINGENT]

        print(f"Suggestive (p < 0.05): {len(suggestive):,}")
        print(f"Stringent (p < 0.001): {len(stringent):,}")

        # Save filtered results
        suggestive.to_csv(OUTPUT_DIR / 'therapy_variants_all.csv', index=False)
        stringent.to_csv(OUTPUT_DIR / 'therapy_variants_stringent.csv', index=False)
        print(f"\nSaved: {OUTPUT_DIR / 'therapy_variants_all.csv'}")
        print(f"Saved: {OUTPUT_DIR / 'therapy_variants_stringent.csv'}")

    # Save cohort summary
    stats_df = pd.DataFrame(cohort_stats)
    stats_df.to_csv(OUTPUT_DIR / 'cohort_summary.csv', index=False)
    print(f"Saved: {OUTPUT_DIR / 'cohort_summary.csv'}")

    # Identify variants appearing in multiple cohorts
    if 'ID' in combined.columns or ('chrom' in combined.columns and 'pos' in combined.columns):
        if 'ID' in combined.columns:
            variant_counts = combined.groupby('ID')['cohort'].nunique()
        else:
            combined['variant_id'] = combined['chrom'].astype(str) + '_' + combined['pos'].astype(str)
            variant_counts = combined.groupby('variant_id')['cohort'].nunique()

        multi_cohort = variant_counts[variant_counts >= 2]
        print(f"\nVariants in ≥2 cohorts: {len(multi_cohort):,}")

        if len(multi_cohort) > 0:
            print("\nTop variants by cohort count:")
            for var, count in multi_cohort.nlargest(10).items():
                print(f"  {var}: {count} cohorts")

    print("\nDone!")


if __name__ == '__main__':
    main()
