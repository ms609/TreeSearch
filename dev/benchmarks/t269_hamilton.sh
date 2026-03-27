#!/bin/bash
#SBATCH --job-name=t269-interleave
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t269_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t269_%j.err

# T-269: Fine-grained sectorial interleaving benchmark
#
# 5 configs × 4 datasets × 5 seeds × {30s, 60s} = 200 runs × ~45s avg ≈ 2.5h
#
# Usage:
#   sbatch t269_hamilton.sh          # 30s budget
#   sbatch t269_hamilton.sh 60       # 60s budget

TIMEOUT=${1:-30}

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t269_results

mkdir -p "$LIB"
mkdir -p "$OUTDIR"
mkdir -p /nobackup/$USER/TreeSearch/logs

echo "=== T-269 Hamilton job ==="
echo "Timeout: ${TIMEOUT}s"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Build and install from latest cpp-search
cd "$REPO" || exit 1
git fetch origin cpp-search
git pull origin cpp-search
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
export R_LIBS_USER="$LIB"
Rscript "$REPO/dev/benchmarks/bench_t269_interleaving.R" "$TIMEOUT" "$OUTDIR"

echo ""
echo "Completed: $(date)"
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t269_*.csv 2>/dev/null
