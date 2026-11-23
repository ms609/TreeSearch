#!/bin/bash
# T-364 battery step 2 on Hamilton: probe P, one array task per arm (1-3).
# Purpose is twofold and both halves matter:
#   (a) GATE -- confirm the frozen constraints actually produce complement-rooted
#       Wagner trees in arm 2, else arm 2 cannot express T-384 and question 2 is
#       unanswerable;
#   (b) DISCRIMINATOR -- prove the three libraries the battery will use really
#       differ.  Three silently identical libs yield a clean, meaningless null.
#SBATCH --job-name=t364-probe
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t364probe_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t364probe_%A_%a.err

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

BASE=/nobackup/$USER/TreeSearch
A=${SLURM_ARRAY_TASK_ID:-1}
export NEOTRANS_DIR=/nobackup/$USER/neotrans/inst/matrices
export CAT_CSV=$BASE-t364arm3/dev/benchmarks/mbank_catalogue.csv
export COMMON=$BASE/t364/t364_common.R
export CONS_DIR=$BASE/t364/cons
export OUT_DIR=$BASE/t364/probe
export TS_LIB=$BASE/t364-lib$A
export R_LIBS="$BASE/t364-lib$A:$BASE/lib"
export ARM=$A
mkdir -p "$OUT_DIR"
cd "$BASE/t364" || exit 1
Rscript "$BASE/t364/t364_probe.R"
