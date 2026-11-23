#!/bin/bash
#SBATCH --job-name=t289c-pr-s2
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=6G
#SBATCH --time=3:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t289c_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t289c_%j.err

# T-289c: Prune-reinsert Stage 2 — mbank_X30754 only, 60s budget
#
# Stage 1 (5 datasets × 13 configs × 5 seeds × 30s) verdict:
#   - ≤88t: PR net-negative (replicate cost >> score gain). Not tested here.
#   - 180t: Real signal. pr_c5_d10 most consistent (5/5 seeds, mean −6.6 steps).
#
# Stage 2 grid: 9 configs × 1 dataset × 10 seeds = 90 runs × ~65s ≈ 98 min.
#   SBATCH --time=3:00:00 provides comfortable margin.
#
# Usage:
#   sbatch t289c_stage2_hamilton.sh [timeout_s]
#   Default timeout: 60s

TIMEOUT=${1:-60}

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t289c_results
export R_LIBS="$LIB:${R_LIBS}"

mkdir -p "$LIB" "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== T-289c Hamilton job (PR Stage 2) ==="
echo "Timeout: ${TIMEOUT}s"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Install CRAN dependencies if missing
echo "Checking CRAN dependencies..."
Rscript --no-save -e "
  lib <- '$LIB'
  .libPaths(c(lib, .libPaths()))
  pkgs <- c('abind', 'ape', 'cli', 'colorspace', 'fastmatch', 'Rdpack', 'TreeDist', 'TreeTools')
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need) > 0) {
    message('Installing: ', paste(need, collapse = ', '))
    install.packages(need, lib = lib, repos = 'https://cloud.r-project.org',
                     dependencies = NA, quiet = TRUE)
  } else {
    message('All dependencies present.')
  }
"

# Build and install TreeSearch from cpp-search
cd "$REPO" || exit 1
git fetch origin cpp-search
git pull --ff-only origin cpp-search || git reset --hard origin/cpp-search
echo "Git HEAD: $(git log --oneline -1)"

INSTALLED_VER=$(Rscript --no-save -e \
  ".libPaths(c('$LIB', .libPaths())); cat(as.character(packageVersion('TreeSearch')))" \
  2>/dev/null || echo "none")
REPO_VER=$(grep '^Version:' DESCRIPTION | awk '{print $2}')
echo "Installed: $INSTALLED_VER  Repo: $REPO_VER"

if [ "$INSTALLED_VER" != "$REPO_VER" ]; then
  echo "Rebuilding TreeSearch..."
  rm -f src/*.o src/*.so
  R CMD build --no-build-vignettes --no-manual --no-resave-data .
  R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
  rc=$?
  rm -f TreeSearch_*.tar.gz
  echo "Install exit code: $rc"
  if [ $rc -ne 0 ]; then
    echo "FATAL: install failed"
    exit 1
  fi
else
  echo "TreeSearch already up to date; skipping rebuild."
fi

# Run benchmark
cd "$OUTDIR"
Rscript "$REPO/dev/benchmarks/bench_pr_stage2_mbank.R" "$TIMEOUT" "$OUTDIR"

echo ""
echo "Completed: $(date)"
echo "Results in: $OUTDIR"
ls -lh "$OUTDIR"/t289c_*.csv 2>/dev/null
