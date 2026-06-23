#!/bin/bash
# Broad byte-identity gate (full bundled NA roster) for the incremental
# exact_verify rescore -- merge-prep before default-on. One cell per dataset;
# each runs EW+IW x several starts, TS_NA_INCR off vs on, asserts identical score.
# Array bound = (#datasets - 1); inapplicable.phyData has 30, so 0-29.
#SBATCH --job-name=ts-naincr-ident
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --array=0-29%30
#SBATCH --output=/nobackup/%u/TreeSearch/logs/naincrid_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/naincrid_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-naiw
REPO=/nobackup/$USER/TreeSearch-t29
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB

cd "$REPO" || exit 1
Rscript dev/benchmarks/bench_na_incr_ident_cell.R "$SLURM_ARRAY_TASK_ID"
