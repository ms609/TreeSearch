#!/bin/bash
#SBATCH --job-name=t253probe
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --time=8:00:00
#SBATCH -o /nobackup/pjjg18/reach/logs/t253probe_%j.out
#SBATCH -e /nobackup/pjjg18/reach/logs/t253probe_%j.err
#
# Decides annotate-vs-retract on the t253 gap analysis.  See the .R for the design.
#
# NOT an array, deliberately: the three steps are SEQUENTIAL (prep feeds both
# addition arms, and both arms feed score), and the whole job is ~150 bare
# `AdditionTree` calls -- an array would need a barrier for no gain.
#
# The archived March engine is the point of this job.  It needs r/4.5.1 (4.4.1
# refuses it: "built under R version 4.5.1") and its own dep path for Rcpp/TreeTools.
set -u
module load r/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

OUT_DIR=/nobackup/pjjg18/reach/t253probe
PROBE=/nobackup/pjjg18/reach/t253_wagner_era_probe.R

CURLIB="/nobackup/pjjg18/curlib:/nobackup/pjjg18/TreeSearch/lib"
# lib-t252 = TreeSearch 2.0.0 packaged 2026-03-27 10:12, i.e. 32 min before
# t252_mbank_30s_20260327_1044.csv was written.  lib-baseline supplies TreeTools 2.2.0
# (and Rcpp, which curlib also needs from TreeSearch/lib).
T252LIB="/nobackup/pjjg18/TreeSearch/lib-t252:/nobackup/pjjg18/ts-bench/lib-baseline"

mkdir -p "$OUT_DIR" /nobackup/pjjg18/reach/logs

echo "=== STEP 1/4: prep (one preprocessing, shared by both arms) ==="
R_LIBS="$CURLIB" \
STEP=prep OUT_DIR="$OUT_DIR" \
NEOTRANS_DIR=/nobackup/pjjg18/neotrans/inst/matrices \
CAT_CSV=/nobackup/pjjg18/reach/mbank_catalogue.csv \
  Rscript "$PROBE" || { echo "PREP FAILED"; exit 1; }

echo "=== STEP 2/4: AdditionTree under the MARCH engine (lib-t252) ==="
R_LIBS="$T252LIB" \
STEP=addition ENGINE=t252 OUT_DIR="$OUT_DIR" N_SEEDS=3 \
  Rscript "$PROBE"

echo "=== STEP 3/4: AdditionTree under the CURRENT engine (curlib) ==="
R_LIBS="$CURLIB" \
STEP=addition ENGINE=cur OUT_DIR="$OUT_DIR" N_SEEDS=3 \
  Rscript "$PROBE"

echo "=== STEP 4/4: score EVERY tree with ONE scorer (curlib) ==="
R_LIBS="$CURLIB" \
STEP=score OUT_DIR="$OUT_DIR" \
  Rscript "$PROBE"

echo "=== done; results in $OUT_DIR ==="
