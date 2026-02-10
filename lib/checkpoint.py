"""
Checkpoint management for resumable variant processing.
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from .log import log


class CheckpointManager:
    """Manage checkpoints for resumable processing with prediction storage."""

    def __init__(self, checkpoint_path: Path, predictions_path: Optional[Path] = None):
        self.checkpoint_path = checkpoint_path
        if predictions_path is not None:
            self.predictions_path = predictions_path
        else:
            stem = checkpoint_path.stem
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

    def mark_processed(self, variant_id: str, index: int, prediction: Optional[Dict] = None):
        """Mark variant as processed and optionally store prediction."""
        self.data['processed_variants'].append(variant_id)
        self.data['last_index'] = index
        if prediction is not None:
            self.predictions[variant_id] = prediction

    def add_error(self, variant_id: str, error: str):
        """Record an error."""
        self.data['errors'].append({
            'variant_id': variant_id,
            'error': error,
            'timestamp': datetime.now().isoformat(),
        })
