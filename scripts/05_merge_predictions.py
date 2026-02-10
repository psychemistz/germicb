#!/usr/bin/env python3
"""
Merge Chunked AlphaGenome Predictions
======================================
After running 02_query_alphagenome.py in parallel chunks via SLURM array,
this script merges per-chunk prediction JSONs into a single unified h5ad.

Input:
- data/therapy_variants_all.csv (or _stringent.csv)
- results/alphagenome/predictions_chunk{0..N}.json

Output:
- results/alphagenome/alphagenome_predictions.h5ad (unified AnnData)

Usage:
  python scripts/05_merge_predictions.py
  python scripts/05_merge_predictions.py --stringent
  python scripts/05_merge_predictions.py --num-chunks 10
"""

import json
import argparse
from pathlib import Path

import pandas as pd

# Reuse functions from the query script
from importlib.util import spec_from_file_location, module_from_spec

# ==============================================================================
# Configuration
# ==============================================================================

PROJECT_DIR = Path('/data/parks34/projects/4germicb')
DATA_DIR = PROJECT_DIR / 'data'
RESULTS_DIR = PROJECT_DIR / 'results' / 'alphagenome'
SCRIPTS_DIR = PROJECT_DIR / 'scripts'


def load_query_module():
    """Import save_predictions_h5ad and prepare_variants from 02_query_alphagenome.py."""
    spec = spec_from_file_location(
        'query_alphagenome',
        SCRIPTS_DIR / '02_query_alphagenome.py',
    )
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    parser = argparse.ArgumentParser(description='Merge chunked AlphaGenome predictions')
    parser.add_argument('--stringent', action='store_true',
                       help='Use stringent variants (p < 0.001)')
    parser.add_argument('--num-chunks', type=int, default=None,
                       help='Expected number of chunks (auto-detected if not set)')
    args = parser.parse_args()

    query_mod = load_query_module()
    log = query_mod.log

    log("=" * 60)
    log("MERGE CHUNKED ALPHAGENOME PREDICTIONS")
    log("=" * 60)

    # Load full variant set (same order as the chunked runs used)
    if args.stringent:
        input_csv = DATA_DIR / 'therapy_variants_stringent.csv'
    else:
        input_csv = DATA_DIR / 'therapy_variants_all.csv'

    log(f"\nLoading variants: {input_csv}")
    variants_df = pd.read_csv(input_csv)
    variants_df = query_mod.prepare_variants(variants_df)
    log(f"  Total variants: {len(variants_df):,}")

    # Discover chunk files
    if args.num_chunks is not None:
        chunk_files = [
            RESULTS_DIR / f'predictions_chunk{i}.json'
            for i in range(args.num_chunks)
        ]
    else:
        chunk_files = sorted(RESULTS_DIR.glob('predictions_chunk*.json'),
                             key=lambda p: int(p.stem.split('chunk')[1]))

    log(f"\nFound {len(chunk_files)} chunk files:")
    for cf in chunk_files:
        log(f"  {cf.name}")

    if len(chunk_files) == 0:
        log("ERROR: No chunk prediction files found in results/alphagenome/")
        log("  Expected files like: predictions_chunk0.json, predictions_chunk1.json, ...")
        return

    # Merge predictions
    merged = {}
    for cf in chunk_files:
        if not cf.exists():
            log(f"  WARNING: Missing chunk file: {cf}")
            continue
        with open(cf, 'r') as f:
            chunk_preds = json.load(f)
        log(f"  {cf.name}: {len(chunk_preds):,} predictions")
        merged.update(chunk_preds)

    log(f"\nMerged total: {len(merged):,} predictions")

    # Check for missing variants
    all_variant_ids = set(variants_df['variant_id'].tolist())
    predicted_ids = set(merged.keys())
    missing = all_variant_ids - predicted_ids
    extra = predicted_ids - all_variant_ids

    if missing:
        log(f"  Missing predictions for {len(missing)} variants")
    if extra:
        log(f"  Extra predictions not in variant list: {len(extra)}")

    # Collect track info
    all_tracks = set()
    for pred in merged.values():
        if 'tracks' in pred:
            all_tracks.update(pred['tracks'].keys())
    log(f"  Unique tracks across all chunks: {len(all_tracks)}")

    # Save unified h5ad
    output_path = RESULTS_DIR / 'alphagenome_predictions.h5ad'
    query_mod.save_predictions_h5ad(variants_df, merged, output_path)

    # Summary
    log("\n" + "=" * 60)
    log("MERGE COMPLETE")
    log("=" * 60)
    log(f"  Total predictions: {len(merged):,}")
    log(f"  Missing variants:  {len(missing)}")
    log(f"  Track count:       {len(all_tracks)}")
    log(f"  Output:            {output_path}")

    log("\nNext steps:")
    log("  1. Run 03_validate_eqtl.py to validate against DICE/OneK1K")
    log("  2. Run 04_prioritize_variants.py to compute priority scores")


if __name__ == '__main__':
    main()
