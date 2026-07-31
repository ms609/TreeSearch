#!/bin/bash
#SBATCH --job-name=reach-ab
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-125
#SBATCH -o /nobackup/pjjg18/reach/logs/reachab_%A_%a.out
#SBATCH -e /nobackup/pjjg18/reach/logs/reachab_%A_%a.err
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# curlib = TreeSearch@cpp-search (49c7fcf, GATE-FREE) + TreeTools@head; deps from TreeSearch/lib.
# The deltas are passed as dots, so no gated build is needed.
export R_LIBS="/nobackup/pjjg18/curlib:/nobackup/pjjg18/TreeSearch/lib"
TS_LIB="/nobackup/pjjg18/curlib" \
NEOTRANS_DIR="/nobackup/pjjg18/neotrans/inst/matrices" \
CAT_CSV="/nobackup/pjjg18/reach/mbank_catalogue.csv" \
OUT_DIR="/nobackup/pjjg18/reach/ab" \
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}" \
N_SEEDS=5 \
  Rscript /nobackup/pjjg18/reach/reach_escalation_ab.R
