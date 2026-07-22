#!/bin/bash
#SBATCH --job-name=t289d-pr-s3
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=6G
#SBATCH --time=3:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t289d_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t289d_%j.err

# T-289d: Prune-reinsert Stage 3 — new drop criteria (MISSING, COMBINED)
#
# Requires TreeSearch >= commit 1ce5e12e (feat: MISSING+COMBINED criteria).
#
# Grid: 8 configs × 1 dataset × 10 seeds × 60s ≈ 87 min.
# SBATCH --time=3:00:00 provides comfortable margin.
#
# Usage:
#   sbatch t289d_stage3_hamilton.sh [timeout_s]
#   Default: 60s

TIMEOUT=${1:-60}

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t289d_results
export R_LIBS="$LIB:${R_LIBS}"

mkdir -p "$LIB" "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== T-289d Hamilton job (PR Stage 3 — new criteria) ==="
echo "Timeout: ${TIMEOUT}s"
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

# Always rebuild — Stage 3 requires the new MISSING/COMBINED criteria
# (commit 1ce5e12e on cpp-search).
cd "$REPO" || exit 1
git fetch origin cpp-search
git pull --ff-only origin cpp-search || git reset --hard origin/cpp-search
echo "Git HEAD: $(git log --oneline -1)"

echo "Rebuilding TreeSearch (new criteria require recompile)..."
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
Rscript "$REPO/dev/benchmarks/bench_pr_stage3_mbank.R" "$TIMEOUT" "$OUTDIR"

echo ""
echo "Completed: $(date)"
ls -lh "$OUTDIR"/t289d_*.csv 2>/dev/null
