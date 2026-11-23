#!/bin/bash
# Hard-tail effort panel: one task per (dataset x seed), all six arms per node.
# Consumes the $LIB built by hamilton_na_hardtail_build.sh.
# 6 matrices x 10 seeds = 60 tasks (0-59).
#
# NO git operations here -- 150 concurrent tasks racing on one index killed 114
# of them last time, and `reset --hard` mutates the checkout surviving tasks are
# reading.  The build job pins the repo; the `Git HEAD:` line below is the
# evidence, and must be identical across all cells before they are pooled.
#SBATCH --job-name=ts-nahard
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
# Worst cell is Zanol2014: effort +2 is rung 5 (1000 replicates at `thorough`
# provisioning, hit target doubled) with certification ON, ~14 h on its own.
# Still well inside `shared`'s 72 h, so no `long` partition.
#SBATCH --time=48:00:00
#SBATCH --array=0-59
#SBATCH --output=/nobackup/%u/TreeSearch/logs/nahard_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/nahard_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-nahard
REPO=/nobackup/$USER/TreeSearch-nahard
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/na_hardtail_partials
export TS_SEEDS="1 2 3 4 5 6 7 8 9 10"
export TS_CONCAVITY=Inf
export TS_DATASETS="Zanol2014 Aguado2009 Zhu2013 Wortley2006 Geisler2001 Aria2015"

mkdir -p "$PARTIAL_DIR" /nobackup/$USER/TreeSearch/logs
cd "$REPO" || exit 1
echo "Git HEAD: $(git log --oneline -1)"
Rscript dev/benchmarks/bench_na_hardtail_cell.R "$SLURM_ARRAY_TASK_ID"
