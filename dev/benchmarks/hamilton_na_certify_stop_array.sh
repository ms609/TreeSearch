#!/bin/bash
# NA-certification under the SHIPPED STOPPING RULES: one task per (dataset x seed),
# all six arms on the same node.  Reuses the $LIB built by
# hamilton_na_certify_build.sh -- no rebuild needed, the mechanism is unchanged.
#
# Panel 1 (hamilton_na_certify_array.sh) neutralised targetHits/stopPatience to
# isolate the budget, and found the gate loses at 8 replicates and wins at matched
# wall.  Production is neither: maxSeconds = 0, maxReplicates = 96,
# targetHits = max(10, ntax/5).  This panel leaves every stopping rule as shipped,
# which is the only regime that can justify changing a preset default.
#
# 30 datasets x 5 seeds = 150 tasks (0-149).
#
# No %N concurrency cap: let SLURM and fairshare pace the array (no-array-throttle).
#SBATCH --job-name=ts-nastop
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
# Headroom for the worst cell: A_hits3 on Zanol2014 runs certification for up to
# the full 96-replicate cap at 3x the default hit target.  Well under `shared`'s
# 72 h, so no `long` partition (hamilton-partition-choice).
#SBATCH --time=24:00:00
#SBATCH --array=0-149
#SBATCH --output=/nobackup/%u/TreeSearch/logs/nastop_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/nastop_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-nacert
REPO=/nobackup/$USER/TreeSearch-nacert
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/na_certify_stop_partials
export TS_SEEDS="1 2 3 4 5"
export TS_CONCAVITY=Inf
# TS_DATASETS unset => all 30 bundled inapplicable matrices.

mkdir -p "$PARTIAL_DIR" /nobackup/$USER/TreeSearch/logs
cd "$REPO" || exit 1
# NO git operations here.  The first submission ran `git fetch && git reset --hard`
# inside the array task, so 150 concurrent tasks raced on one index and 114 died
# on `.git/index.lock`.  Worse than the failure: `reset --hard` mutates the
# checkout that the surviving tasks are reading, so a "successful" cell could have
# run a half-rewritten script.  Sync the repo ONCE from the login node (or the
# build job) before submitting, and keep it PINNED for the life of the panel so
# every cell provably runs the same code -- the `Git HEAD:` line below is the
# evidence, and it must be identical across all cells before they are pooled.
echo "Git HEAD: $(git log --oneline -1)"
Rscript dev/benchmarks/bench_na_certify_stop_cell.R "$SLURM_ARRAY_TASK_ID"
