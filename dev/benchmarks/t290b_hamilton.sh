#!/bin/bash
#SBATCH --job-name=t290b-phase
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t290b_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t290b_%j.err

# T-290b: Phase profiling across scoring modes and weightings
# 6 datasets x 2 scoring x 2 weightings x 2 strategies x 3 seeds = 144 runs
# @ 30s each = ~72 min + build/overhead
# Estimated total: ~90 min

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t290_results

mkdir -p "$OUTDIR"
mkdir -p /nobackup/$USER/TreeSearch/logs

echo "=== T-290b Phase Profiling ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

cd "$REPO" || exit 1
echo "Git HEAD: $(git log --oneline -1)"
echo ""

# Verify bench scripts and library exist
for f in bench_datasets.R t290b_phase_profile.R; do
  if [ ! -f "$REPO/dev/benchmarks/$f" ]; then
    echo "FATAL: $f not found"
    exit 1
  fi
done

Rscript -e ".libPaths(c('$LIB', .libPaths())); library(TreeSearch); message('TreeSearch ', packageVersion('TreeSearch'))"
if [ $? -ne 0 ]; then
  echo "FATAL: TreeSearch not loadable from $LIB"
  exit 1
fi

NEOTRANS=/nobackup/$USER/neotrans/inst/matrices
if [ ! -d "$NEOTRANS" ]; then
  echo "FATAL: neotrans matrices not found at $NEOTRANS"
  exit 1
fi
echo "Neotrans matrices: $(ls $NEOTRANS | wc -l) files"
echo ""

cd "$REPO"
export R_LIBS_USER="$LIB"
Rscript dev/benchmarks/t290b_phase_profile.R "$OUTDIR" 2>&1

echo ""
echo "Completed: $(date)"
ls -la "$OUTDIR"/t290b_*.csv 2>/dev/null
