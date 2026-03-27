#!/bin/bash
#SBATCH --job-name=t252-mbank
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t252_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t252_%j.err

# T-252: MorphoBank training-set baseline benchmark
# Phase 1: Install dependencies
# Phase 2: Install TreeSearch
# Phase 3: Run 25 matrices x 3 budgets x 5 seeds = 375 runs

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t252_results

mkdir -p "$LIB" "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== T-252 MorphoBank Training-Set Benchmark ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Phase 1: Install R dependencies
echo "=== Phase 1: Installing R dependencies ==="
export R_LIBS_USER="$LIB"
Rscript -e "
  .libPaths(c('$LIB', .libPaths()))
  needed <- c('Rcpp', 'ape', 'TreeTools', 'TreeDist', 'Rdpack',
              'cli', 'fastmatch', 'abind', 'colorspace')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cat('Installing:', paste(missing, collapse = ', '), '\n')
    install.packages(missing, lib = '$LIB',
                     repos = 'https://cloud.r-project.org', Ncpus = 1)
  } else {
    cat('All dependencies already installed\n')
  }
  # Verify
  ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop('Still missing: ', paste(needed[!ok], collapse = ', '))
  }
  cat('All', length(needed), 'dependencies OK\n')
" 2>&1
rc=$?
if [ $rc -ne 0 ]; then
  echo "FATAL: dependency installation failed"
  exit 1
fi

# Phase 2: Install TreeSearch
echo ""
echo "=== Phase 2: Installing TreeSearch ==="
cd "$REPO" || exit 1
git pull --ff-only origin cpp-search 2>/dev/null || true
echo "Git HEAD: $(git log --oneline -1)"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
echo "Install exit code: $rc"
rm -f TreeSearch_*.tar.gz

if [ $rc -ne 0 ]; then
  echo "FATAL: TreeSearch install failed"
  exit 1
fi

# Verify neotrans
NEOTRANS=/nobackup/$USER/neotrans/inst/matrices
if [ ! -d "$NEOTRANS" ] || [ "$(ls $NEOTRANS | wc -l)" -eq 0 ]; then
  echo "FATAL: neotrans matrices not found or empty at $NEOTRANS"
  exit 1
fi
echo "Neotrans matrices: $(ls $NEOTRANS | wc -l) files"

# Phase 3: Run benchmark
echo ""
echo "=== Phase 3: Running benchmark ==="
cd "$REPO"
Rscript dev/benchmarks/bench_t252_mbank_training.R "$OUTDIR" 2>&1

echo ""
echo "=== Completed: $(date) ==="
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t252_*.csv 2>/dev/null
