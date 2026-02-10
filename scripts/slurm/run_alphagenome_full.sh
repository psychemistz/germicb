#!/bin/bash
#SBATCH --job-name=germicb_ag
#SBATCH --partition=norm
#SBATCH --array=0-9
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --output=/data/parks34/projects/4germicb/logs/alphagenome_chunk%a_%A.out
#SBATCH --error=/data/parks34/projects/4germicb/logs/alphagenome_chunk%a_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=seongyong.park@nih.gov

# ==============================================================================
# Run AlphaGenome Predictions — SLURM Array Job (10 parallel chunks)
# ==============================================================================
#
# Each array task processes ~1/10 of variants using batch predict_variants().
# After all tasks complete, submit the merge + validation job:
#
#   ARRAY_JOB_ID=$(sbatch --parsable scripts/slurm/run_alphagenome_full.sh)
#   sbatch --dependency=afterok:${ARRAY_JOB_ID} scripts/slurm/run_merge_validate.sh
#
# Or to resume failed chunks:
#   sbatch --array=3,7 scripts/slurm/run_alphagenome_full.sh
#
# ==============================================================================

set -e

NUM_CHUNKS=10
CHUNK_INDEX=${SLURM_ARRAY_TASK_ID}

echo "=============================================="
echo "AlphaGenome Batch Predictions — Chunk ${CHUNK_INDEX}/${NUM_CHUNKS}"
echo "=============================================="
echo "Job ID:       ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node:         ${SLURMD_NODENAME}"
echo "Start time:   $(date)"
echo "Time limit:   2 days"
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
    echo "ERROR: ALPHAGENOME_API_KEY not set"
    echo "Set it in your ~/.bashrc or submit with:"
    echo "  ALPHAGENOME_API_KEY=your_key sbatch ..."
    exit 1
fi

# Run AlphaGenome predictions for this chunk
echo "Starting AlphaGenome batch predictions..."
echo "  Chunk: ${CHUNK_INDEX} of ${NUM_CHUNKS}"
echo "  Batch size: 50, Max workers: 5, Seq length: 100kb"
echo ""

python scripts/02_query_alphagenome.py \
  --chunk-index ${CHUNK_INDEX} \
  --num-chunks ${NUM_CHUNKS} \
  --max-workers 5 \
  --batch-size 50 \
  --seq-length 100kb

echo ""
echo "=============================================="
echo "Chunk ${CHUNK_INDEX} completed: $(date)"
echo "=============================================="
