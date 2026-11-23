#!/bin/bash
# Job-array panel: one task per (dataset x seed) cell, consuming the shared $LIB
# from hamilton_build_once.sh (submit with --dependency=afterok:<buildjob>).
# Tune --array upper bound to (n_datasets * n_seeds - 1); %N caps concurrency.
# DISPATCH-UNTESTED template — bench_cell.R is validated locally.
#SBATCH --job-name=ts-panel
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --array=0-29%32
#SBATCH --output=/nobackup/%u/TreeSearch/logs/panel_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/panel_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

LIB=/nobackup/$USER/TreeSearch/lib
REPO=/nobackup/$USER/TreeSearch-a
export R_LIBS_USER=$LIB
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/panel_partials
export TS_REPS=20
export TS_SEEDS="1 2 3 4 5"
export TS_DATASETS="Wortley2006 Eklund2004 Zanol2014 Zhu2013 Giles2015 Dikow2009"

cd "$REPO" || exit 1
Rscript dev/benchmarks/bench_cell.R "$SLURM_ARRAY_TASK_ID"
