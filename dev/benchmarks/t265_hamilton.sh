#!/bin/bash
#SBATCH --job-name=t265-regression
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t265_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t265_%j.err

# T-265: Per-replicate quality regression diagnosis
# 3 configs x 9 datasets x 5 seeds x 120s = ~135 runs x ~120s = ~4.5 hours

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t265_results

mkdir -p "$LIB"
mkdir -p "$OUTDIR"

echo "=== T-265 Hamilton job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Build and install from latest cpp-search
cd "$REPO" || exit 1
echo "Git HEAD: $(git log --oneline -1)"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
echo "Install exit code: $rc"
rm -f TreeSearch_*.tar.gz

if [ $rc -ne 0 ]; then
  echo "FATAL: install failed"
  exit 1
fi

# Run benchmark
cd "$OUTDIR"
Rscript -e ".libPaths(c('$LIB', .libPaths()))" \
  "$REPO/dev/benchmarks/bench_t265_regression.R" 120 "$OUTDIR"

echo ""
echo "Completed: $(date)"
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t265_*.csv 2>/dev/null
