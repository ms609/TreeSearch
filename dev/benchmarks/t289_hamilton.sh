#!/bin/bash
#SBATCH --job-name=t289-prune-ri
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=8:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t289_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t289_%j.err

# T-289: Prune-reinsert perturbation benchmark
#
# Stage 1: 13 configs × 5 datasets × 5 seeds × 30s ≈ 325 runs × ~30s ≈ 2.7h
# Stage 2: ~10 configs × 5 datasets × 5 seeds × {30s,60s} ≈ 500 runs × ~45s ≈ 6.3h
#
# Usage:
#   sbatch t289_hamilton.sh           # runs stage 1 (30s)
#   sbatch t289_hamilton.sh 2 30      # stage 2, 30s budget
#   sbatch t289_hamilton.sh 2 60      # stage 2, 60s budget

STAGE=${1:-1}
TIMEOUT=${2:-30}

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t289_results
export R_LIBS="$LIB:${R_LIBS}"

mkdir -p "$LIB"
mkdir -p "$OUTDIR"
mkdir -p /nobackup/$USER/TreeSearch/logs

echo "=== T-289 Hamilton job ==="
echo "Stage: $STAGE, Timeout: ${TIMEOUT}s"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Install CRAN dependencies into local lib (if missing)
echo "Checking/installing CRAN dependencies..."
Rscript --no-save -e "
  lib <- '$LIB'
  .libPaths(c(lib, .libPaths()))
  pkgs <- c('abind', 'ape', 'cli', 'colorspace', 'fastmatch', 'Rdpack', 'TreeDist', 'TreeTools')
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need) > 0) {
    message('Installing: ', paste(need, collapse = ', '))
    install.packages(need, lib = lib, repos = 'https://cloud.r-project.org',
                     dependencies = NA, quiet = TRUE)  # NA = Imports+Depends only
  } else {
    message('All dependencies present.')
  }
"

# Build and install from latest cpp-search
cd "$REPO" || exit 1
git fetch origin cpp-search
git pull --ff-only origin cpp-search || git reset --hard origin/cpp-search
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
Rscript "$REPO/dev/benchmarks/bench_prune_reinsert.R" "$STAGE" "$TIMEOUT" "$OUTDIR"

echo ""
echo "Completed: $(date)"
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t289_*.csv 2>/dev/null
