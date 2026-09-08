#!/bin/bash
#SBATCH --job-name=2island-sweep
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --array=0-31
#SBATCH --output=/nobackup/%u/TreeSearch/logs/2island_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/2island_%A_%a.err
#
# Array sweep: tune "thorough" for two-island completeness + settle wagnerStarts.
# Each task runs a deterministic stripe of the cell grid (driver shards by
# SLURM_ARRAY_TASK_ID over SLURM_ARRAY_TASK_COUNT). Submit AFTER the build:
#   bid=$(sbatch --parsable hamilton_two_island_build.sh)
#   sbatch --dependency=afterok:$bid hamilton_two_island_sweep.sh
# Then aggregate: Rscript dev/benchmarks/hamilton_two_island_aggregate.R
#
# Grid (defaults): PART A 6 variants x {90,180}s x 30 seeds = 360 cells;
#                  PART B 6 variants x {60,120}s x 15 seeds x 4 datasets = 720.
# ~1080 cells, ~32 CPU-h total; 32 shards -> ~1 CPU-h/task.

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-2island
LIB=/nobackup/$USER/TreeSearch/lib2island
OUTDIR=/nobackup/$USER/TreeSearch/two_island_results
mkdir -p "$OUTDIR"
export R_LIBS="$LIB:${R_LIBS}"

# Guard: build must have completed.
[ -f "$LIB/.build_done" ] || { echo "FATAL: build sentinel missing; run build job first"; exit 2; }

echo "=== shard $SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_COUNT === $(date) | node $(hostname)"

cd "$OUTDIR"
TS_LIB="$LIB" REPO="$REPO" OUTDIR="$OUTDIR" \
  NSEED_A=30 NSEED_B=15 BUDGETS_A="90 180" BUDGETS_B="60 120" \
  SHARD="$SLURM_ARRAY_TASK_ID" NSHARD="$SLURM_ARRAY_TASK_COUNT" \
  Rscript "$REPO/dev/benchmarks/hamilton_two_island_sweep.R"

echo "=== shard $SLURM_ARRAY_TASK_ID done === $(date)"
