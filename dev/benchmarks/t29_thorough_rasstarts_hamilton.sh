#!/bin/bash
#SBATCH --job-name=t29-thorough-ras
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t29_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t29_%j.err

# Task #29: full-search time-matched gate for rasStarts=3 in the AUTO-SELECTED
# `thorough` preset.  Local indicative + rss-only time-matched both favour
# rasStarts=3 (Zanol/Zhu, 5-8 steps); this confirms on authoritative wall-clock
# over the full thorough pipeline before flipping thorough's default (intensive
# already adopted it, commit e69765f3).
#
# Grid: 4 datasets x rasStarts{1,3} x budgets{60,120}s x 10 seeds = 160 runs.
# Expected ~4h; 8h limit gives margin.

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-t29
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t29_results
export R_LIBS="$LIB:${R_LIBS}"
mkdir -p "$LIB" "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== Task #29: thorough rasStarts=1 vs 3 (full-search, time-matched) ==="
echo "Job ID: $SLURM_JOB_ID | Node: $(hostname) | Started: $(date)"

Rscript --no-save -e "
  lib <- '$LIB'; .libPaths(c(lib, .libPaths()))
  pkgs <- c('abind','ape','cli','colorspace','fastmatch','Rdpack','TreeDist','TreeTools')
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) install.packages(need, lib = lib,
        repos = 'https://cloud.r-project.org', dependencies = NA, quiet = TRUE)
"

# Build from origin/cpp-search (carries the rasStarts lever + build_ras_sector fix
# + this driver).  Clone if absent.
if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/ms609/TreeSearch.git "$REPO"
fi
cd "$REPO" || exit 1
git fetch origin cpp-search && git reset --hard origin/cpp-search
echo "Git HEAD: $(git log --oneline -1)"

echo "Rebuilding TreeSearch..."
rm -f src/*.o src/*.so src/*.dll
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz; rc=$?
rm -f TreeSearch_*.tar.gz
if [ $rc -ne 0 ]; then echo "FATAL: install failed"; exit 1; fi

cd "$OUTDIR"
TS_LIB="$LIB" OUTDIR="$OUTDIR" NSEED=10 BUDGETS="60 120" \
  Rscript "$REPO/dev/benchmarks/hamilton_thorough_rasstarts.R"

echo "Completed: $(date)"
ls -lh "$OUTDIR"/thorough_rasstarts.csv 2>/dev/null
