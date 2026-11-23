#!/bin/bash
#SBATCH --job-name=t289b-brazeau
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=8:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t289b_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t289b_%j.err

# T-289b: Prune-reinsert benchmark — Brazeau (default) scoring
#
# Parallel companion to t289_hamilton.sh (Fitch/EW mode).
# Uses TreeSearch's default Brazeau et al. (2019) inapplicable scoring.
# Shares the same build artifact as the Fitch job — no rebuild needed
# if t289_hamilton.sh has already installed TreeSearch in $LIB.
#
# Stage 1: 13 configs x 5 datasets x 5 seeds x 30s ≈ 325 runs ≈ 2.7h
# Stage 2: ~10 configs x 5 datasets x 5 seeds x {30s,60s} ≈ 500 runs ≈ 6.3h
#
# Usage:
#   sbatch t289b_brazeau_hamilton.sh           # stage 1, 30s
#   sbatch t289b_brazeau_hamilton.sh 2 30      # stage 2, 30s
#   sbatch t289b_brazeau_hamilton.sh 2 60      # stage 2, 60s

STAGE=${1:-1}
TIMEOUT=${2:-30}

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t289b_results
export R_LIBS="$LIB:${R_LIBS}"

mkdir -p "$LIB"
mkdir -p "$OUTDIR"
mkdir -p /nobackup/$USER/TreeSearch/logs

echo "=== T-289b Hamilton job (Brazeau scoring) ==="
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
                     dependencies = NA, quiet = TRUE)
  } else {
    message('All dependencies present.')
  }
"

# Build and install from latest cpp-search
# (Skip rebuild if TreeSearch is already installed and up to date;
#  rebuild if the Fitch job hasn't run yet or if HEAD has moved.)
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
  echo "Install exit code: $rc"
  rm -f TreeSearch_*.tar.gz
  if [ $rc -ne 0 ]; then
    echo "FATAL: install failed"
    exit 1
  fi
else
  echo "TreeSearch already up to date; skipping rebuild."
fi

# Run benchmark
cd "$OUTDIR"
Rscript "$REPO/dev/benchmarks/bench_prune_reinsert_brazeau.R" "$STAGE" "$TIMEOUT" "$OUTDIR"

echo ""
echo "Completed: $(date)"
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t289b_*.csv 2>/dev/null
