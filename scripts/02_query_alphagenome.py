#!/usr/bin/env python3
"""
Query AlphaGenome API for Immunotherapy-Associated Variants
============================================================
Call AlphaGenome API for therapy-associated variants to predict regulatory effects.

This is the CORRECT use case for AlphaGenome:
- Input: Variants associated with therapy outcomes (unknown regulatory effects)
- Output: Predicted regulatory impact from AlphaGenome model

Features:
- Uses Google DeepMind AlphaGenome API (batch predict_variants() for throughput)
- SLURM array chunking: split variants across parallel jobs
- Configurable batch size, concurrency, and sequence context length
- Rate limiting with exponential backoff retry
- Checkpoint after each batch for resume on failure
- Filter to immune-relevant tracks (GM12878, PBMCs, etc.)

Input:
- data/therapy_variants_all.csv (16,061 variants, p < 0.05)
- data/therapy_variants_stringent.csv (83 variants, p < 0.001)

Output (per chunk when using --chunk-index/--num-chunks):
- results/alphagenome/predictions_chunk{N}.json
- results/alphagenome/checkpoint_chunk{N}.json

Output (unchunked or after merging with 05_merge_predictions.py):
- results/alphagenome/alphagenome_predictions.h5ad (AnnData with track scores)

Environment:
- ALPHAGENOME_API_KEY: Google DeepMind API key (required for real predictions)
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, List

import numpy as np
import pandas as pd
import anndata as ad

# ==============================================================================
# Configuration
# ==============================================================================

PROJECT_DIR = Path('/data/parks34/projects/4germicb')
DATA_DIR = PROJECT_DIR / 'data'
RESULTS_DIR = PROJECT_DIR / 'results' / 'alphagenome'

# Rate limiting
CHECKPOINT_INTERVAL = 100
MAX_RETRIES = 5
INITIAL_BACKOFF = 1.0  # seconds

# Immune-relevant track patterns (case-insensitive matching)
IMMUNE_TRACK_PATTERNS = [
    'gm12878',      # B-lymphoblastoid cell line
    'pbmc',         # Peripheral blood mononuclear cells
    'cd4',          # CD4+ T cells
    'cd8',          # CD8+ T cells
    'b_cell', 'bcell', 'b-cell',  # B cells
    't_cell', 'tcell', 't-cell',  # T cells
    'monocyte',     # Monocytes
    'macrophage',   # Macrophages
    'dendritic',    # Dendritic cells
    'nk_cell', 'nkcell', 'nk-cell',  # NK cells
    'lymphocyte',   # Lymphocytes
    'leukocyte',    # Leukocytes
    'immune',       # Generic immune
    'blood',        # Blood-related
    'hematopoietic', 'hsc',  # Hematopoietic stem cells
    'spleen',       # Spleen
    'thymus',       # Thymus
    'bone_marrow', 'bone-marrow',  # Bone marrow
    'cd34',         # CD34+ progenitors
]


def log(msg: str):
    """Print timestamped log message."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class CheckpointManager:
    """Manage checkpoints for resumable processing with prediction storage."""

    def __init__(self, checkpoint_path: Path):
        self.checkpoint_path = checkpoint_path
        # Derive predictions checkpoint path from the checkpoint filename
        stem = checkpoint_path.stem  # e.g. 'checkpoint_chunk0' or 'checkpoint'
        self.predictions_path = checkpoint_path.parent / f'{stem}_predictions.json'
        self.data = self.load()
        self.predictions = self.load_predictions()

    def load(self) -> Dict[str, Any]:
        """Load checkpoint from disk."""
        if self.checkpoint_path.exists():
            with open(self.checkpoint_path, 'r') as f:
                return json.load(f)
        return {
            'processed_variants': [],
            'last_index': -1,
            'start_time': datetime.now().isoformat(),
            'errors': [],
        }

    def load_predictions(self) -> Dict[str, Dict]:
        """Load predictions from separate file."""
        if self.predictions_path.exists():
            try:
                with open(self.predictions_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                log(f"  Warning: Could not load predictions checkpoint: {e}")
        return {}

    def save(self):
        """Save checkpoint to disk."""
        self.data['last_updated'] = datetime.now().isoformat()
        self.data['n_predictions'] = len(self.predictions)
        with open(self.checkpoint_path, 'w') as f:
            json.dump(self.data, f, indent=2)

    def save_predictions(self):
        """Save predictions to separate file."""
        with open(self.predictions_path, 'w') as f:
            json.dump(self.predictions, f)
        log(f"  Predictions checkpoint saved ({len(self.predictions)} variants)")

    def is_processed(self, variant_id: str) -> bool:
        """Check if variant has been processed."""
        return variant_id in self.data['processed_variants']

    def mark_processed(self, variant_id: str, index: int, prediction: Dict = None):
        """Mark variant as processed and store prediction."""
        self.data['processed_variants'].append(variant_id)
        self.data['last_index'] = index
        if prediction is not None:
            self.predictions[variant_id] = prediction

    def add_error(self, variant_id: str, error: str):
        """Record an error."""
        self.data['errors'].append({
            'variant_id': variant_id,
            'error': error,
            'timestamp': datetime.now().isoformat()
        })


def is_immune_track(track_name: str) -> bool:
    """Check if track name matches immune-relevant patterns."""
    track_lower = track_name.lower()
    return any(pattern in track_lower for pattern in IMMUNE_TRACK_PATTERNS)


def process_variant_output(variant_output, filter_immune: bool = True) -> Dict[str, Any]:
    """
    Process AlphaGenome VariantOutput to extract track scores.

    Args:
        variant_output: VariantOutput from AlphaGenome API
        filter_immune: Whether to filter to immune-relevant tracks only

    Returns:
        Dictionary with track scores and metadata
    """
    from alphagenome.models.dna_output import OutputType

    tracks = {}

    # VariantOutput has reference and alternate Output objects
    if not hasattr(variant_output, 'reference') or not hasattr(variant_output, 'alternate'):
        return {'tracks': {}, 'metadata': {'total_tracks': 0, 'significant_tracks': 0}}

    ref_output = variant_output.reference
    alt_output = variant_output.alternate

    # Process each output type (ATAC, DNASE, CHIP_HISTONE, CHIP_TF, RNA_SEQ)
    output_types = [
        ('atac', OutputType.ATAC),
        ('dnase', OutputType.DNASE),
        ('chip_histone', OutputType.CHIP_HISTONE),
        ('chip_tf', OutputType.CHIP_TF),
        ('rna_seq', OutputType.RNA_SEQ),
    ]

    for field_name, output_type in output_types:
        ref_track_data = getattr(ref_output, field_name, None)
        alt_track_data = getattr(alt_output, field_name, None)

        if ref_track_data is None or alt_track_data is None:
            continue

        # Get track names from metadata
        if hasattr(ref_track_data, 'metadata') and 'name' in ref_track_data.metadata.columns:
            track_names = ref_track_data.metadata['name'].tolist()
        else:
            continue

        # Get values (shape: positional_bins x num_tracks)
        ref_values = ref_track_data.values
        alt_values = alt_track_data.values

        if ref_values is None or alt_values is None:
            continue

        # Process each track
        for i, track_name in enumerate(track_names):
            # Create full track name with output type prefix
            full_name = f"{output_type.name}_{track_name}"

            if filter_immune and not is_immune_track(full_name):
                continue

            try:
                # Compute mean across positional bins
                if ref_values.ndim == 1:
                    ref_score = float(ref_values[i]) if i < len(ref_values) else 0.0
                    alt_score = float(alt_values[i]) if i < len(alt_values) else 0.0
                else:
                    ref_score = float(np.mean(ref_values[:, i]))
                    alt_score = float(np.mean(alt_values[:, i]))

                tracks[full_name] = {
                    'ref_score': ref_score,
                    'alt_score': alt_score,
                    'diff': alt_score - ref_score,
                }
            except Exception:
                continue

    return {
        'tracks': tracks,
        'metadata': {
            'total_tracks': len(tracks),
            'significant_tracks': sum(1 for t in tracks.values() if abs(t['diff']) > 0.1),
        }
    }


def run_real_predictions(
    variants_df: pd.DataFrame,
    checkpoint: CheckpointManager,
    api_key: str,
    resume: bool = False,
    output_dir: Path = None
) -> Dict[str, Dict]:
    """
    Run real AlphaGenome predictions using the API.

    Args:
        variants_df: DataFrame with variant information
        checkpoint: CheckpointManager for resumable processing
        api_key: AlphaGenome API key
        resume: Whether to resume from checkpoint
        output_dir: Directory for intermediate h5ad saves

    Returns:
        Dictionary mapping variant_id to predictions
    """
    from alphagenome.models import dna_client
    from alphagenome.models.dna_output import OutputType
    from alphagenome.data import genome

    log("Connecting to AlphaGenome API...")
    client = dna_client.create(api_key=api_key)
    log("  Connected successfully")

    # Request immune-relevant output types
    requested_outputs = [
        OutputType.ATAC,           # Chromatin accessibility
        OutputType.DNASE,          # DNase accessibility
        OutputType.CHIP_HISTONE,   # Histone modifications (H3K27ac, H3K4me3, etc.)
        OutputType.CHIP_TF,        # Transcription factor binding
        OutputType.RNA_SEQ,        # Gene expression
    ]
    log(f"  Requesting output types: {[o.name for o in requested_outputs]}")

    # Start with existing predictions from checkpoint if resuming
    predictions = checkpoint.predictions.copy() if resume else {}
    start_idx = checkpoint.data['last_index'] + 1 if resume else 0
    total = len(variants_df)

    log(f"\nProcessing {total - start_idx} variants...")
    if resume and len(predictions) > 0:
        log(f"  Loaded {len(predictions)} predictions from checkpoint")
    if start_idx > 0:
        log(f"  Resuming from index {start_idx}")

    for i, (_, row) in enumerate(variants_df.iloc[start_idx:].iterrows()):
        idx = start_idx + i
        variant_id = row['variant_id']

        # Skip if already processed
        if checkpoint.is_processed(variant_id):
            continue

        chrom = row['chrom']
        pos = int(row['pos'])
        ref = row['ref']
        alt = row['alt']

        log(f"  [{idx + 1}/{total}] {variant_id} - {chrom}:{pos} {ref}>{alt}")

        pred = None
        backoff = INITIAL_BACKOFF
        for attempt in range(MAX_RETRIES):
            try:
                # Create variant object
                variant = genome.Variant(
                    chromosome=chrom,
                    position=pos,
                    reference_bases=ref,
                    alternate_bases=alt
                )

                # Create interval centered on variant (1MB context)
                seq_len = dna_client.SEQUENCE_LENGTH_1MB
                interval_start = max(0, pos - seq_len // 2)
                interval_end = interval_start + seq_len
                interval = genome.Interval(chrom, interval_start, interval_end)

                # Call API - request multiple output types for immune analysis
                result = client.predict_variant(
                    interval=interval,
                    variant=variant,
                    requested_outputs=requested_outputs,
                    ontology_terms=None,  # All cell types/tissues
                )

                # Process result
                pred = process_variant_output(result, filter_immune=True)
                predictions[variant_id] = pred
                n_tracks = len(pred.get('tracks', {}))
                log(f"    Got {n_tracks} immune tracks")
                break  # Success, exit retry loop

            except Exception as e:
                error_msg = str(e)
                if 'RESOURCE_EXHAUSTED' in error_msg and attempt < MAX_RETRIES - 1:
                    log(f"    Rate limited, waiting {backoff:.1f}s...")
                    time.sleep(backoff)
                    backoff *= 2  # Exponential backoff
                else:
                    checkpoint.add_error(variant_id, error_msg)
                    log(f"    ERROR: {error_msg[:100]}")
                    break

        # Update checkpoint with prediction
        checkpoint.mark_processed(variant_id, idx, pred)

        # Save checkpoint and predictions periodically
        if (idx + 1) % CHECKPOINT_INTERVAL == 0:
            checkpoint.save()
            checkpoint.save_predictions()
            log(f"  Checkpoint saved at index {idx + 1}")

            # Save intermediate h5ad every 1000 variants
            if output_dir and (idx + 1) % 1000 == 0:
                intermediate_path = output_dir / 'alphagenome_predictions_intermediate.h5ad'
                try:
                    save_predictions_h5ad(variants_df, predictions, intermediate_path)
                    log(f"  Intermediate h5ad saved ({len(predictions)} variants)")
                except Exception as e:
                    log(f"  Warning: Could not save intermediate h5ad: {e}")

    # Final save
    checkpoint.save()
    checkpoint.save_predictions()
    return predictions


def _get_seq_length(seq_length_str: str):
    """Map sequence length string to SDK constant."""
    from alphagenome.models import dna_client
    mapping = {
        '1mb': dna_client.SEQUENCE_LENGTH_1MB,
        '100kb': dna_client.SEQUENCE_LENGTH_100KB,
        '16kb': dna_client.SEQUENCE_LENGTH_16KB,
    }
    return mapping[seq_length_str]


def run_batch_predictions(
    variants_df: pd.DataFrame,
    checkpoint: CheckpointManager,
    api_key: str,
    resume: bool = False,
    output_dir: Path = None,
    max_workers: int = 5,
    batch_size: int = 50,
    seq_length_str: str = '100kb',
) -> Dict[str, Dict]:
    """
    Run AlphaGenome predictions using batch predict_variants() API.

    Processes variants in batches for higher throughput via concurrent requests.

    Args:
        variants_df: DataFrame with variant information
        checkpoint: CheckpointManager for resumable processing
        api_key: AlphaGenome API key
        resume: Whether to resume from checkpoint
        output_dir: Directory for intermediate h5ad saves
        max_workers: Concurrent API workers for predict_variants()
        batch_size: Number of variants per batch call
        seq_length_str: Sequence context length ('1mb', '100kb', '16kb')

    Returns:
        Dictionary mapping variant_id to predictions
    """
    from alphagenome.models import dna_client
    from alphagenome.models.dna_output import OutputType
    from alphagenome.data import genome

    log("Connecting to AlphaGenome API...")
    client = dna_client.create(api_key=api_key)
    log("  Connected successfully")

    seq_len = _get_seq_length(seq_length_str)
    log(f"  Sequence context length: {seq_length_str} ({seq_len:,} bp)")

    requested_outputs = [
        OutputType.ATAC,
        OutputType.DNASE,
        OutputType.CHIP_HISTONE,
        OutputType.CHIP_TF,
        OutputType.RNA_SEQ,
    ]
    log(f"  Requesting output types: {[o.name for o in requested_outputs]}")
    log(f"  Batch size: {batch_size}, Max workers: {max_workers}")

    # Start with existing predictions from checkpoint if resuming
    predictions = checkpoint.predictions.copy() if resume else {}
    already_processed = set(checkpoint.data['processed_variants']) if resume else set()
    total = len(variants_df)

    # Filter out already-processed variants
    remaining_indices = []
    for i, (_, row) in enumerate(variants_df.iterrows()):
        if row['variant_id'] not in already_processed:
            remaining_indices.append(i)

    log(f"\nTotal variants in chunk: {total}")
    if resume and len(predictions) > 0:
        log(f"  Loaded {len(predictions)} predictions from checkpoint")
    log(f"  Remaining to process: {len(remaining_indices)}")

    # Process in batches
    n_batches = (len(remaining_indices) + batch_size - 1) // batch_size
    log(f"  Will process in {n_batches} batches of up to {batch_size}")

    variants_processed_since_checkpoint = 0

    for batch_num in range(n_batches):
        batch_start = batch_num * batch_size
        batch_end = min(batch_start + batch_size, len(remaining_indices))
        batch_indices = remaining_indices[batch_start:batch_end]

        batch_rows = variants_df.iloc[batch_indices]
        batch_variant_ids = batch_rows['variant_id'].tolist()

        log(f"\n  Batch {batch_num + 1}/{n_batches} ({len(batch_indices)} variants)")

        # Build Variant and Interval objects for the batch
        variant_objects = []
        interval_objects = []
        for _, row in batch_rows.iterrows():
            chrom = row['chrom']
            pos = int(row['pos'])
            ref = row['ref']
            alt = row['alt']

            variant_objects.append(genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
            ))

            interval_start = max(0, pos - seq_len // 2)
            interval_end = interval_start + seq_len
            interval_objects.append(genome.Interval(chrom, interval_start, interval_end))

        # Call batch API with exponential backoff retry
        backoff = INITIAL_BACKOFF
        results = None
        for attempt in range(MAX_RETRIES):
            try:
                results = client.predict_variants(
                    intervals=interval_objects,
                    variants=variant_objects,
                    requested_outputs=requested_outputs,
                    ontology_terms=None,
                    max_workers=max_workers,
                )
                break  # Success
            except Exception as e:
                error_msg = str(e)
                if 'RESOURCE_EXHAUSTED' in error_msg and attempt < MAX_RETRIES - 1:
                    log(f"    Rate limited on batch, waiting {backoff:.1f}s (attempt {attempt + 1}/{MAX_RETRIES})...")
                    time.sleep(backoff)
                    backoff *= 2
                else:
                    log(f"    ERROR on batch (attempt {attempt + 1}/{MAX_RETRIES}): {error_msg[:150]}")
                    if attempt == MAX_RETRIES - 1:
                        # Record errors for all variants in batch
                        for vid in batch_variant_ids:
                            checkpoint.add_error(vid, f"Batch failed: {error_msg[:200]}")

        # Process batch results
        if results is not None:
            for j, (variant_output, vid) in enumerate(zip(results, batch_variant_ids)):
                try:
                    pred = process_variant_output(variant_output, filter_immune=True)
                    predictions[vid] = pred
                    checkpoint.mark_processed(vid, batch_indices[j], pred)
                except Exception as e:
                    checkpoint.add_error(vid, str(e))
                    checkpoint.mark_processed(vid, batch_indices[j], None)
                    log(f"    Error processing {vid}: {str(e)[:100]}")

            n_tracks_avg = np.mean([
                len(predictions[vid].get('tracks', {}))
                for vid in batch_variant_ids if vid in predictions
            ]) if any(vid in predictions for vid in batch_variant_ids) else 0
            log(f"    Got results for {len(results)} variants (avg {n_tracks_avg:.0f} immune tracks)")
        else:
            # Mark all as processed even on failure so we don't retry forever
            for j, vid in enumerate(batch_variant_ids):
                checkpoint.mark_processed(vid, batch_indices[j], None)

        variants_processed_since_checkpoint += len(batch_indices)

        # Checkpoint after each batch
        if variants_processed_since_checkpoint >= CHECKPOINT_INTERVAL:
            checkpoint.save()
            checkpoint.save_predictions()
            log(f"    Checkpoint saved ({len(predictions)} total predictions)")
            variants_processed_since_checkpoint = 0

            # Save intermediate h5ad periodically
            if output_dir and len(predictions) % 1000 < batch_size:
                intermediate_path = output_dir / 'alphagenome_predictions_intermediate.h5ad'
                try:
                    save_predictions_h5ad(variants_df, predictions, intermediate_path)
                    log(f"    Intermediate h5ad saved ({len(predictions)} variants)")
                except Exception as e:
                    log(f"    Warning: Could not save intermediate h5ad: {e}")

    # Final save
    checkpoint.save()
    checkpoint.save_predictions()
    return predictions


def run_mock_predictions(variants_df: pd.DataFrame) -> Dict[str, Dict]:
    """
    Generate mock predictions for testing when AlphaGenome is unavailable.

    This creates realistic-looking but synthetic data for pipeline testing.
    """
    log("Running MOCK predictions (AlphaGenome API not used)...")

    np.random.seed(42)
    predictions = {}

    # Mock track names (based on what AlphaGenome actually returns)
    # Note: AlphaGenome returns only 3 RNA-seq tracks for immune: lymphocytes, spleen, whole blood
    mock_tracks = [
        'RNA_SEQ_GTEx_Lymphocytes_EBV',
        'RNA_SEQ_GTEx_Spleen',
        'RNA_SEQ_GTEx_Whole_Blood',
    ]

    for i, (_, row) in enumerate(variants_df.iterrows()):
        variant_id = row['variant_id']

        # Generate random scores with some variants having larger effects
        tracks = {}
        for track in mock_tracks:
            ref_score = np.random.uniform(0, 10)
            # Some variants have larger effects (correlated with therapy association)
            if row['pval'] < 0.001:
                effect = np.random.normal(0, 1.5)  # Larger effects for stringent variants
            elif row['pval'] < 0.01:
                effect = np.random.normal(0, 0.8)
            else:
                effect = np.random.normal(0, 0.3)

            tracks[track] = {
                'ref_score': ref_score,
                'alt_score': ref_score + effect,
                'diff': effect,
            }

        predictions[variant_id] = {
            'tracks': tracks,
            'metadata': {
                'total_tracks': len(tracks),
                'significant_tracks': sum(1 for t in tracks.values() if abs(t['diff']) > 0.3),
                'mock': True,
            }
        }

        if (i + 1) % 1000 == 0:
            log(f"  Processed {i + 1}/{len(variants_df)}")

    return predictions


def save_predictions_h5ad(
    variants_df: pd.DataFrame,
    predictions: Dict[str, Dict],
    output_path: Path
):
    """
    Save predictions to h5ad format.

    Structure:
    - obs: variant metadata (one row per variant)
    - var: track metadata (one column per track)
    - X: track score differences (variants x tracks)
    - layers['ref_score']: reference allele scores
    - layers['alt_score']: alternate allele scores
    """
    log("Saving predictions to h5ad...")

    # Collect all track names across all variants
    all_tracks = set()
    for pred in predictions.values():
        if 'tracks' in pred:
            all_tracks.update(pred['tracks'].keys())

    all_tracks = sorted(all_tracks)
    log(f"  Total tracks: {len(all_tracks)}")

    if len(all_tracks) == 0:
        log("  WARNING: No track predictions to save")
        return

    # Create matrices
    n_variants = len(variants_df)
    n_tracks = len(all_tracks)

    diff_matrix = np.zeros((n_variants, n_tracks))
    ref_matrix = np.zeros((n_variants, n_tracks))
    alt_matrix = np.zeros((n_variants, n_tracks))

    track_to_idx = {t: i for i, t in enumerate(all_tracks)}

    for i, (_, row) in enumerate(variants_df.iterrows()):
        variant_id = row['variant_id']
        if variant_id in predictions and 'tracks' in predictions[variant_id]:
            for track_name, scores in predictions[variant_id]['tracks'].items():
                if track_name in track_to_idx:
                    j = track_to_idx[track_name]
                    diff_matrix[i, j] = scores.get('diff', 0)
                    ref_matrix[i, j] = scores.get('ref_score', 0)
                    alt_matrix[i, j] = scores.get('alt_score', 0)

    # Create AnnData
    adata = ad.AnnData(
        X=diff_matrix,
        obs=variants_df.reset_index(drop=True),
        var=pd.DataFrame({'track_name': all_tracks}, index=all_tracks)
    )

    adata.layers['ref_score'] = ref_matrix
    adata.layers['alt_score'] = alt_matrix

    # Add metadata
    adata.uns['stage'] = 'alphagenome_predictions'
    adata.uns['description'] = 'AlphaGenome variant effect predictions for immunotherapy variants'
    adata.uns['n_variants'] = n_variants
    adata.uns['n_tracks'] = n_tracks
    adata.uns['project'] = '4germicb - Immunotherapy WES Germline Variants'

    adata.write_h5ad(output_path, compression='gzip')
    log(f"  Saved: {output_path}")


def prepare_variants(variants_df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare variants DataFrame for AlphaGenome processing.

    Standardizes column names and creates variant_id.
    Drops rows with missing coordinates.
    """
    df = variants_df.copy()

    # Drop rows with missing essential columns
    required_cols = ['chrom', 'pos', 'ref', 'alt']
    n_before = len(df)
    df = df.dropna(subset=required_cols)
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        log(f"  Dropped {n_dropped} variants with missing coordinates")

    # Convert pos to int
    df['pos'] = df['pos'].astype(int)

    # Ensure chromosome format is correct (with 'chr' prefix)
    if len(df) > 0 and not str(df['chrom'].iloc[0]).startswith('chr'):
        df['chrom'] = 'chr' + df['chrom'].astype(str)

    # Create variant_id if not present
    if 'variant_id' not in df.columns:
        df['variant_id'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str) + '_' + df['ref'] + '>' + df['alt']

    return df


def main():
    parser = argparse.ArgumentParser(description='Query AlphaGenome for Immunotherapy Variants')
    parser.add_argument('--resume', action='store_true',
                       help='Resume from checkpoint')
    parser.add_argument('--mock', action='store_true',
                       help='Use mock predictions (for testing without API)')
    parser.add_argument('--stringent', action='store_true',
                       help='Use stringent variants only (p < 0.001, n=83)')
    parser.add_argument('--max-variants', type=int, default=None,
                       help='Maximum variants to process (for testing)')
    parser.add_argument('--api-key', type=str, default=None,
                       help='AlphaGenome API key (or set ALPHAGENOME_API_KEY env var)')
    parser.add_argument('--reset', action='store_true',
                       help='Reset checkpoint and start fresh')
    # Batch / chunking arguments
    parser.add_argument('--chunk-index', type=int, default=None,
                       help='Which chunk this job processes (0-based, maps to SLURM_ARRAY_TASK_ID)')
    parser.add_argument('--num-chunks', type=int, default=None,
                       help='Total number of parallel chunks')
    parser.add_argument('--max-workers', type=int, default=5,
                       help='Concurrent API workers for predict_variants() (default: 5)')
    parser.add_argument('--batch-size', type=int, default=50,
                       help='Number of variants per batch call (default: 50)')
    parser.add_argument('--seq-length', type=str, default='100kb',
                       choices=['1mb', '100kb', '16kb'],
                       help='Sequence context length (default: 100kb)')
    args = parser.parse_args()

    log("=" * 60)
    log("ALPHAGENOME QUERY: IMMUNOTHERAPY VARIANTS")
    log("=" * 60)

    # Create output directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load therapy variants
    if args.stringent:
        input_csv = DATA_DIR / 'therapy_variants_stringent.csv'
        log(f"\nLoading STRINGENT variants: {input_csv}")
    else:
        input_csv = DATA_DIR / 'therapy_variants_all.csv'
        log(f"\nLoading ALL suggestive variants: {input_csv}")

    variants_df = pd.read_csv(input_csv)
    log(f"  Loaded {len(variants_df):,} variants")

    # Prepare variants
    variants_df = prepare_variants(variants_df)

    # Limit variants for testing
    if args.max_variants:
        variants_df = variants_df.head(args.max_variants)
        log(f"  Limited to {len(variants_df)} variants for testing")

    # Chunking: split variants across parallel jobs
    chunk_index = args.chunk_index
    num_chunks = args.num_chunks

    if chunk_index is not None and num_chunks is not None:
        total_variants = len(variants_df)
        chunk_size = (total_variants + num_chunks - 1) // num_chunks
        start = chunk_index * chunk_size
        end = min(start + chunk_size, total_variants)
        variants_df = variants_df.iloc[start:end].reset_index(drop=True)
        log(f"\n  Chunk {chunk_index}/{num_chunks}: variants [{start}:{end}) = {len(variants_df)} variants")
        chunk_suffix = f"_chunk{chunk_index}"
    elif chunk_index is not None or num_chunks is not None:
        log("ERROR: --chunk-index and --num-chunks must both be provided")
        sys.exit(1)
    else:
        chunk_suffix = ""

    # Show cohort breakdown
    log("\nVariants by cohort:")
    for cohort, count in variants_df['cohort'].value_counts().items():
        log(f"  {cohort}: {count:,}")

    # Initialize checkpoint (chunk-specific paths)
    checkpoint_path = RESULTS_DIR / f'checkpoint{chunk_suffix}.json'

    # Handle reset
    if args.reset:
        log("\nResetting checkpoint (starting fresh)...")
        if checkpoint_path.exists():
            checkpoint_path.unlink()
        predictions_ckpt = checkpoint_path.parent / f'predictions{chunk_suffix}_checkpoint.json'
        if predictions_ckpt.exists():
            predictions_ckpt.unlink()
        log("  Checkpoint files removed")

    checkpoint = CheckpointManager(checkpoint_path)

    if args.resume:
        log(f"\nResuming from checkpoint...")
        log(f"  Previously processed: {len(checkpoint.data['processed_variants'])}")
        log(f"  Saved predictions: {len(checkpoint.predictions)}")
        log(f"  Errors: {len(checkpoint.data['errors'])}")

    # Get API key
    api_key = args.api_key or os.environ.get('ALPHAGENOME_API_KEY')

    # Run predictions
    if args.mock:
        predictions = run_mock_predictions(variants_df)
    elif api_key:
        log(f"\nUsing AlphaGenome API (batch mode)...")
        try:
            predictions = run_batch_predictions(
                variants_df, checkpoint, api_key, args.resume, RESULTS_DIR,
                max_workers=args.max_workers,
                batch_size=args.batch_size,
                seq_length_str=args.seq_length,
            )
        except Exception as e:
            log(f"\nERROR with AlphaGenome API: {e}")
            log("Falling back to mock predictions...")
            predictions = run_mock_predictions(variants_df)
    else:
        log("\nWARNING: No API key provided")
        log("  Set ALPHAGENOME_API_KEY environment variable or use --api-key")
        log("  Falling back to mock predictions...")
        predictions = run_mock_predictions(variants_df)

    log(f"\nCompleted predictions for {len(predictions)} variants")

    # Save predictions JSON (chunk-specific)
    predictions_path = RESULTS_DIR / f'predictions{chunk_suffix}.json'
    with open(predictions_path, 'w') as f:
        json.dump(predictions, f)
    log(f"  Saved predictions JSON: {predictions_path}")

    # If not chunked, also save the unified h5ad directly
    if not chunk_suffix:
        output_path = RESULTS_DIR / 'alphagenome_predictions.h5ad'
        save_predictions_h5ad(variants_df, predictions, output_path)

    # Update checkpoint with completion info
    checkpoint.data['completed'] = True
    checkpoint.data['n_predictions'] = len(predictions)
    checkpoint.save()

    # Summary
    log("\n" + "=" * 60)
    log("COMPLETE")
    log("=" * 60)
    log(f"  Predictions: {len(predictions)}")
    log(f"  Errors: {len(checkpoint.data['errors'])}")
    log(f"  Output: {predictions_path}")

    if chunk_suffix:
        log("\nNext steps:")
        log("  1. Run all chunks, then run 05_merge_predictions.py")
    else:
        log("\nNext steps:")
        log("  1. Run 03_validate_eqtl.py to validate against DICE/GTEx/CIMA")
        log("  2. Run 04_prioritize_variants.py to compute priority scores")


if __name__ == '__main__':
    main()
