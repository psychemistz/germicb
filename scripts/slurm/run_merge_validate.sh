#!/bin/bash
#SBATCH --job-name=germicb_merge
#SBATCH --partition=norm
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --output=/data/parks34/projects/4germicb/logs/merge_validate_%j.out
#SBATCH --error=/data/parks34/projects/4germicb/logs/merge_validate_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=seongyong.park@nih.gov

# ==============================================================================
# Merge Chunked Predictions + Run Downstream Validation & Prioritization
# ==============================================================================
#
# Run after all array tasks from run_alphagenome_full.sh complete:
#
#   ARRAY_JOB_ID=$(sbatch --parsable scripts/slurm/run_alphagenome_full.sh)
#   sbatch --dependency=afterok:${ARRAY_JOB_ID} scripts/slurm/run_merge_validate.sh
#
# ==============================================================================

set -e

echo "=============================================="
echo "Merge & Validate AlphaGenome Predictions"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $SLURMD_NODENAME"
echo "Start:  $(date)"
echo ""

# Setup environment
source ~/bin/myconda
conda activate secactpy

cd /data/parks34/projects/4germicb

# Step 1: Merge chunked predictions into unified h5ad
echo "Step 1: Merging chunked predictions..."
python scripts/05_merge_predictions.py --num-chunks 10

# Step 2: Validate against eQTL databases
echo ""
echo "Step 2: Running eQTL validation..."
python scripts/03_validate_eqtl.py

# Step 3: Prioritize variants
echo ""
echo "Step 3: Running variant prioritization..."
python scripts/04_prioritize_variants.py

echo ""
echo "=============================================="
echo "All steps completed: $(date)"
echo "=============================================="
