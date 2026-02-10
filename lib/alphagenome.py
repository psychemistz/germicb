"""
AlphaGenome output processing and immune track filtering.
"""

from typing import Any, Dict

import numpy as np

# Immune-relevant track patterns (case-insensitive matching)
IMMUNE_TRACK_PATTERNS = [
    'gm12878',      # B-lymphoblastoid cell line
    'pbmc',         # Peripheral blood mononuclear cells
    'cd4',          # CD4+ T cells
    'cd8',          # CD8+ T cells
    'b_cell', 'bcell', 'b-cell',
    't_cell', 'tcell', 't-cell',
    'monocyte',
    'macrophage',
    'dendritic',
    'nk_cell', 'nkcell', 'nk-cell',
    'lymphocyte',
    'leukocyte',
    'immune',
    'blood',
    'hematopoietic', 'hsc',
    'spleen',
    'thymus',
    'bone_marrow', 'bone-marrow',
    'cd34',
]


def is_immune_track(track_name: str) -> bool:
    """Check if track name matches immune-relevant patterns."""
    track_lower = track_name.lower()
    return any(pattern in track_lower for pattern in IMMUNE_TRACK_PATTERNS)


def filter_immune_tracks(tracks: Dict[str, Any]) -> Dict[str, Any]:
    """Filter a tracks dictionary to only immune-relevant entries."""
    return {k: v for k, v in tracks.items() if is_immune_track(k)}


def process_variant_output(variant_output, filter_immune: bool = True) -> Dict[str, Any]:
    """
    Process AlphaGenome VariantOutput to extract track scores.

    Args:
        variant_output: VariantOutput from AlphaGenome API
        filter_immune: Whether to filter to immune-relevant tracks only

    Returns:
        Dictionary with 'tracks' and 'metadata' keys
    """
    from alphagenome.models.dna_output import OutputType

    tracks = {}

    if not hasattr(variant_output, 'reference') or not hasattr(variant_output, 'alternate'):
        return {'tracks': {}, 'metadata': {'total_tracks': 0, 'significant_tracks': 0}}

    ref_output = variant_output.reference
    alt_output = variant_output.alternate

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

        if hasattr(ref_track_data, 'metadata') and 'name' in ref_track_data.metadata.columns:
            track_names = ref_track_data.metadata['name'].tolist()
        else:
            continue

        ref_values = ref_track_data.values
        alt_values = alt_track_data.values

        if ref_values is None or alt_values is None:
            continue

        for i, track_name in enumerate(track_names):
            full_name = f"{output_type.name}_{track_name}"

            if filter_immune and not is_immune_track(full_name):
                continue

            try:
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
