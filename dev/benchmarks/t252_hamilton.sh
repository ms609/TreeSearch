#!/bin/bash
#SBATCH --job-name=t252-mbank
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t252_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t252_%j.err

# T-252: MorphoBank training-set baseline benchmark
# 25 matrices x 3 budgets (30/60/120s) x 5 seeds = 375 runs
# Estimated: ~5 hours

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t252_results

mkdir -p "$LIB"
mkdir -p "$OUTDIR"
mkdir -p /nobackup/$USER/TreeSearch/logs

echo "=== T-252 MorphoBank Training-Set Benchmark ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Build and install from latest cpp-search
cd "$REPO" || exit 1
git pull --ff-only origin cpp-search 2>/dev/null || true
echo "Git HEAD: $(git log --oneline -1)"
echo ""

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

# Verify neotrans corpus is available
NEOTRANS=/nobackup/$USER/neotrans/inst/matrices
if [ ! -d "$NEOTRANS" ]; then
  echo "FATAL: neotrans matrices not found at $NEOTRANS"
  echo "Clone with: cd /nobackup/$USER && git clone <neotrans-repo>"
  exit 1
fi
echo "Neotrans matrices: $(ls $NEOTRANS | wc -l) files"
echo ""

# Run benchmark
cd "$REPO"
export R_LIBS_USER="$LIB"
Rscript dev/benchmarks/bench_t252_mbank_training.R "$OUTDIR" 2>&1

echo ""
echo "Completed: $(date)"
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t252_*.csv 2>/dev/null
