#!/bin/bash
# IW mission-wall A/B job-array: one task per (dataset x seed) cell, each running
# BOTH arms (opts off vs on) so the wall ratio controls for node variance.
# Consumes the shared $LIB from hamilton_build_once.sh (the build must include the
# IW opts + the kill-switches TS_IW_NOX4/TS_IW_NODIRTY — i.e. the branch carrying
# these changes, NOT plain cpp-search until merged).
# Submit: sbatch --dependency=afterok:<buildjob> hamilton_iw_wall_array.sh
# Array upper bound = n_datasets * n_seeds - 1 (here 4*5-1=19); %N caps concurrency.
# Merge partials with hamilton_merge.sh (or just cat the iwcell_*.csv).
#SBATCH --job-name=ts-iw-wall
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --array=0-19%20
#SBATCH --output=/nobackup/%u/TreeSearch/logs/iwwall_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/iwwall_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

LIB=${TS_LIB:-/nobackup/$USER/TreeSearch/lib}
REPO=${TS_REPO:-/nobackup/$USER/TreeSearch-t29}
export R_LIBS_USER=$LIB
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/iw_wall_partials
export TS_REPS=10
export TS_SEEDS="1 2 3 4 5"
# Pure-IW: recoded inside the cell runner ("-"->"?"). Mission roster, IW-relevant.
export TS_DATASETS="Wortley2006 Zanol2014 Giles2015 Dikow2009"

cd "$REPO" || exit 1
Rscript dev/benchmarks/bench_iw_wall_cell.R "$SLURM_ARRAY_TASK_ID"
