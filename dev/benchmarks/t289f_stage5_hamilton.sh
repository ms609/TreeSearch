#!/bin/bash
#SBATCH --job-name=t289f-pr-nni
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t289f_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t289f_%j.err

# T-289f: Prune-reinsert Stage 5 — NNI full-tree polish cost reduction
#
# Compares: baseline (no PR) vs pr_nni (NNI polish) vs pr_tbr (TBR polish)
# on the same 5 large-tree datasets as Stage 4 (131-206 tips).
#
# Builds from feature/tbr-batch (contains pruneReinsertNni parameter).
#
# Grid: 5 datasets x 3 configs x 2 budgets x 10 seeds = 300 runs
# Expected wall time: ~4-6h; 8h limit provides comfortable margin.

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t289f_results
export R_LIBS="$LIB:${R_LIBS}"

mkdir -p "$LIB" "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== T-289f Hamilton job (PR Stage 5 — NNI Polish) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Install CRAN dependencies if missing
Rscript --no-save -e "
  lib <- '$LIB'
  .libPaths(c(lib, .libPaths()))
  pkgs <- c('abind', 'ape', 'cli', 'colorspace', 'fastmatch', 'Rdpack', 'TreeDist', 'TreeTools')
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need) > 0) {
    message('Installing: ', paste(need, collapse = ', '))
    install.packages(need, lib = lib, repos = 'https://cloud.r-project.org',
                     dependencies = NA, quiet = TRUE)
  } else { message('All dependencies present.') }
"

# Rebuild from cpp-search (pruneReinsertNni merged via PR #238)
cd "$REPO" || exit 1
git fetch origin cpp-search
git reset --hard origin/cpp-search
echo "Git HEAD: $(git log --oneline -1)"

echo "Rebuilding TreeSearch..."
rm -f src/*.o src/*.so src/*.dll
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz
echo "Install exit code: $rc"
if [ $rc -ne 0 ]; then
  echo "FATAL: install failed"
  exit 1
fi

# Run benchmark
cd "$OUTDIR"
Rscript "$REPO/dev/benchmarks/bench_pr_stage5_nni.R" "$OUTDIR"

echo ""
echo "Completed: $(date)"
ls -lh "$OUTDIR"/t289f_*.csv 2>/dev/null
