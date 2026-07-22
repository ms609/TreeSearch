#!/bin/bash
# NA-IW x4 mission-wall A/B job-array: one task per (dataset x seed) cell, each
# running BOTH arms (x4 off vs on) so the wall ratio controls for node variance.
# Consumes the $LIB built by hamilton_iw_na_build.sh (the NA-IW-x4 branch).
# Submit: sbatch --dependency=afterok:<buildjob> hamilton_iw_na_array.sh
# Array upper bound = n_datasets * n_seeds - 1 (here 4*5-1=19); %N caps concurrency.
# Native NA datasets (NOT recoded) + XPIWE (MaximizeParsimony default) = production.
#SBATCH --job-name=ts-naiw-wall
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=3:00:00
#SBATCH --array=0-19%20
#SBATCH --output=/nobackup/%u/TreeSearch/logs/naiwwall_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/naiwwall_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-naiw
REPO=/nobackup/$USER/TreeSearch-t29
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/iw_na_wall_partials
export TS_REPS=8
export TS_SEEDS="1 2 3 4 5"
# Native NA, largest trees spanning the NA-fraction range (7-64%).
export TS_DATASETS="Dikow2009 Giles2015 Zanol2014 Zhu2013"

cd "$REPO" || exit 1
Rscript dev/benchmarks/bench_iw_na_wall_cell.R "$SLURM_ARRAY_TASK_ID"
