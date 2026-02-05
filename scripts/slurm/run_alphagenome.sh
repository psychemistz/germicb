#!/bin/bash
#SBATCH --job-name=germicb_alphagenome
#SBATCH --partition=norm
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --output=/data/parks34/projects/4germicb/logs/alphagenome_%j.out
#SBATCH --error=/data/parks34/projects/4germicb/logs/alphagenome_%j.err

# ==============================================================================
# Run AlphaGenome Predictions for Immunotherapy Variants
# ==============================================================================
#
# Usage:
#   sbatch scripts/slurm/run_alphagenome.sh [--stringent] [--resume]
#
# Options:
#   --stringent  Use stringent variants only (p < 0.001, n=83)
#   --resume     Resume from checkpoint
#
# Note: Set ALPHAGENOME_API_KEY environment variable before running
# ==============================================================================

set -e

# Parse arguments
EXTRA_ARGS=""
for arg in "$@"; do
    EXTRA_ARGS="$EXTRA_ARGS $arg"
done

echo "=============================================="
echo "AlphaGenome Predictions for Immunotherapy"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Date: $(date)"
echo "Arguments: $EXTRA_ARGS"
echo ""

# Setup environment
source ~/bin/myconda
conda activate secactpy

# Create log directory
mkdir -p /data/parks34/projects/4germicb/logs

# Change to project directory
cd /data/parks34/projects/4germicb

# Check for API key
if [ -z "$ALPHAGENOME_API_KEY" ]; then
    echo "WARNING: ALPHAGENOME_API_KEY not set"
    echo "Running with mock predictions for testing"
    EXTRA_ARGS="$EXTRA_ARGS --mock"
fi

# Run AlphaGenome predictions
echo ""
echo "Running AlphaGenome predictions..."
python scripts/02_query_alphagenome.py $EXTRA_ARGS

echo ""
echo "=============================================="
echo "Completed: $(date)"
echo "=============================================="
