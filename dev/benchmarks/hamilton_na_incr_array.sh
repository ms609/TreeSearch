#!/bin/bash
# NA incremental-rescore mission-wall A/B (paired TS_NA_INCR off/on).
# Needs lib-naiw REBUILT at the commit carrying the incremental (50046bfa) --
# re-run hamilton_iw_na_build.sh first (chain on afterok).
#SBATCH --job-name=ts-naincr-wall
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --array=0-19%20
#SBATCH --output=/nobackup/%u/TreeSearch/logs/naincr_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/naincr_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-naiw
REPO=/nobackup/$USER/TreeSearch-t29
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/na_incr_partials
export TS_REPS=5
export TS_SEEDS="1 2 3 4 5"
export TS_DATASETS="Dikow2009 Giles2015 Zanol2014 Zhu2013"

cd "$REPO" || exit 1
Rscript dev/benchmarks/bench_na_incr_wall_cell.R "$SLURM_ARRAY_TASK_ID"
