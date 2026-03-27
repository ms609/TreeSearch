#!/bin/bash
#SBATCH --job-name=t252-mbank
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t252_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t252_%j.err

# T-252: MorphoBank training-set baseline benchmark (v2 — fixed lib paths)
# 25 matrices x 3 budgets (30/60/120s) x 5 seeds = 375 runs (~5 hours)
#
# Uses ts-bench/lib-baseline for all deps (TreeDist, TreeTools, etc.),
# installs only the fresh TreeSearch build into TreeSearch/lib-t252.

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
FRESH_LIB=/nobackup/$USER/TreeSearch/lib-t252
DEP_LIB=/nobackup/$USER/ts-bench/lib-baseline
OUTDIR=/nobackup/$USER/TreeSearch/t252_results

mkdir -p "$FRESH_LIB" "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== T-252 MorphoBank Training-Set Benchmark v2 ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Fresh lib: $FRESH_LIB"
echo "Dep lib:   $DEP_LIB"
echo ""

# Phase 1: Build and install TreeSearch (deps resolved from DEP_LIB)
echo "=== Building TreeSearch from cpp-search ==="
cd "$REPO" || exit 1
git pull --ff-only origin cpp-search 2>/dev/null || true
echo "Git HEAD: $(git log --oneline -1)"

rm -f src/*.o src/*.so
TMPBUILD=$(mktemp -d)
(cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$REPO")

# Install using both libs so R can find TreeSearch's Imports during install
export R_LIBS="$FRESH_LIB:$DEP_LIB"
R CMD INSTALL --library="$FRESH_LIB" "$TMPBUILD"/TreeSearch_*.tar.gz
rc=$?
rm -rf "$TMPBUILD"
echo "Install exit code: $rc"

if [ $rc -ne 0 ]; then
  echo "FATAL: TreeSearch install failed"
  exit 1
fi

# Verify the install loaded correctly
Rscript -e "
  .libPaths(c('$FRESH_LIB', '$DEP_LIB', .libPaths()))
  library(TreeSearch)
  cat('TreeSearch version:', as.character(packageVersion('TreeSearch')), '\n')
"
rc=$?
if [ $rc -ne 0 ]; then
  echo "FATAL: TreeSearch failed to load"
  exit 1
fi

# Phase 2: Verify neotrans corpus
NEOTRANS=/nobackup/$USER/neotrans/inst/matrices
if [ ! -d "$NEOTRANS" ] || [ "$(ls $NEOTRANS | wc -l)" -eq 0 ]; then
  echo "FATAL: neotrans matrices not found at $NEOTRANS"
  exit 1
fi
echo "Neotrans matrices: $(ls $NEOTRANS | wc -l) files"

# Phase 3: Run benchmark
echo ""
echo "=== Running benchmark ==="
cd "$REPO"
export R_LIBS="$FRESH_LIB:$DEP_LIB"
Rscript dev/benchmarks/bench_t252_mbank_training.R "$OUTDIR" 2>&1

echo ""
echo "=== Completed: $(date) ==="
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t252_*.csv 2>/dev/null
