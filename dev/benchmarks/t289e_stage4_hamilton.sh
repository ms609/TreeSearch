#!/bin/bash
#SBATCH --job-name=t289e-pr-s4
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t289e_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t289e_%j.err

# T-289e: Prune-reinsert Stage 4 — multi-dataset validation
#
# Validates that PR (c=5, d=5%, MISSING) benefit generalises across 5 large-tree
# matrices (131-206 tips) and persists at 120s budget.
#
# Grid: 5 datasets × 2 configs × 2 budgets × 10 seeds = 200 runs
# Expected wall time: ~5h; 8h limit provides comfortable margin.

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t289e_results
export R_LIBS="$LIB:${R_LIBS}"

mkdir -p "$LIB" "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== T-289e Hamilton job (PR Stage 4 — multi-dataset validation) ==="
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

# Rebuild — Stage 4 runs the large preset with PR (commit in cpp-search)
cd "$REPO" || exit 1
git fetch origin cpp-search
git pull --ff-only origin cpp-search || git reset --hard origin/cpp-search
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
Rscript "$REPO/dev/benchmarks/bench_pr_stage4_validation.R" "$OUTDIR"

echo ""
echo "Completed: $(date)"
ls -lh "$OUTDIR"/t289e_*.csv 2>/dev/null
