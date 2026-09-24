#!/bin/bash
# NA-IW x4 ELEMENT-level A/B (dilution-free direct ts_tbr_search climbs).
# Reuses the lib-naiw built by hamilton_iw_na_build.sh (package src unchanged since
# b34df50a -- only benchmark scripts were added), so NO rebuild dependency needed;
# just ensure TreeSearch-t29 is pulled to the commit carrying this runner.
# Submit: sbatch hamilton_iw_na_element_array.sh
#SBATCH --job-name=ts-naiw-elem
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --array=0-19%20
#SBATCH --output=/nobackup/%u/TreeSearch/logs/naiwelem_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/naiwelem_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-naiw
REPO=/nobackup/$USER/TreeSearch-t29
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/iw_na_elem_partials
export TS_REPS=5
export TS_SEEDS="1 2 3 4 5"
export TS_DATASETS="Dikow2009 Giles2015 Zanol2014 Zhu2013"

cd "$REPO" || exit 1
Rscript dev/benchmarks/bench_iw_na_element_cell.R "$SLURM_ARRAY_TASK_ID"
