#!/bin/bash
# DSS size-axis probe panel: one array task per (dataset x seed) manifest row (1..12),
# consuming the dedicated $LIB from hamilton_dss_build.sh. Submit with
# --dependency=afterok:<buildjob>. Confirms the local pilot at n=15 (3 T0 seeds x 5 reps)
# on the CORRECTED axis (wall + reach-to-target; steps/Mcand reported but flagged biased).
#SBATCH --job-name=ts-dss-panel
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=0:40:00
#SBATCH --array=1-12
#SBATCH --output=/nobackup/%u/TreeSearch/logs/dsspanel_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/dsspanel_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

LIB=/nobackup/$USER/TreeSearch/dsslib
DEPLIB=/nobackup/$USER/TreeSearch/lib      # deps (TreeTools/ape/Rcpp/...); my TreeSearch is in $LIB
REPO=/nobackup/$USER/TreeSearch-dss
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB                          # dss_size_probe.R pins TreeSearch via lib.loc=TS_LIB
export T0_DIR=$REPO/dev/benchmarks/missiongate/t0
export OUT_DIR=/nobackup/$USER/TreeSearch/dss_out
export TS_REPS=5
mkdir -p "$OUT_DIR"

cd "$REPO" || exit 1
Rscript dev/benchmarks/missiongate/dss_size_probe.R
