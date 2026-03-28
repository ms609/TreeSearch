#!/bin/bash
#SBATCH --job-name=t290-brazeau
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t290_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t290_%j.err

# T-290: Brazeau-track baseline benchmark
# 20 matrices x 2 weightings (EW, IW k=10) x 2 strategies x 5 seeds x 30s
# = 400 runs @ 30s each = ~3.5 hours + overhead
# Estimated total: ~5 hours

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-a
LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/t290_results

mkdir -p "$LIB"
mkdir -p "$OUTDIR"
mkdir -p /nobackup/$USER/TreeSearch/logs

echo "=== T-290 Brazeau-Track Baseline Benchmark ==="
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

# The bench_brazeau_baseline.R script lives in the TS-TNT-bench worktree.
# Copy it to the repo so it can source bench_datasets.R normally.
# (bench_datasets.R expects to be sourced from the repo root.)
BENCH_SCRIPT="$REPO/dev/benchmarks/bench_brazeau_baseline.R"
if [ ! -f "$BENCH_SCRIPT" ]; then
  # Try fetching from the tnt-bench branch
  cd "$REPO"
  git show feature/tnt-bench:dev/benchmarks/bench_brazeau_baseline.R > "$BENCH_SCRIPT" 2>/dev/null
  git show feature/tnt-bench:dev/benchmarks/bench_datasets.R > "$REPO/dev/benchmarks/bench_datasets.R.tnt" 2>/dev/null
  if [ -f "$REPO/dev/benchmarks/bench_datasets.R.tnt" ]; then
    cp "$REPO/dev/benchmarks/bench_datasets.R.tnt" "$REPO/dev/benchmarks/bench_datasets.R"
    rm "$REPO/dev/benchmarks/bench_datasets.R.tnt"
  fi
  git show feature/tnt-bench:dev/benchmarks/t290_run.R > "$REPO/dev/benchmarks/t290_run.R" 2>/dev/null
fi

if [ ! -f "$BENCH_SCRIPT" ]; then
  echo "FATAL: bench_brazeau_baseline.R not found"
  exit 1
fi

if [ ! -f "$REPO/dev/benchmarks/t290_run.R" ]; then
  echo "FATAL: t290_run.R not found"
  exit 1
fi

# Run benchmark
cd "$REPO"
export R_LIBS_USER="$LIB"
Rscript dev/benchmarks/t290_run.R "$OUTDIR" 2>&1

echo ""
echo "Completed: $(date)"
echo "Results in: $OUTDIR"
ls -la "$OUTDIR"/t290_*.csv 2>/dev/null
