#!/bin/bash
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 01:30:00
#SBATCH --job-name=tsmarch
#SBATCH --output=/nobackup/pjjg18/tsmarch/logs/tsmarch_%j.out
#SBATCH --error=/nobackup/pjjg18/tsmarch/logs/tsmarch_%j.err
# REUSABLE deployment build: TreeSearch (cpp-search) with -march=native for THIS
# machine class (EPYC 7702 / Zen2). Removes the per-block AVX2-call boundary the
# shipped -msse2 build pays -> ~1-10% end-to-end (larger on scoring-bound char-rich
# data). Machine-specific: point search jobs' R_LIBS_USER / lib.loc at $DEPLOY.
# NOT for CRAN (-march=native breaks non-AVX2 CPUs; release keeps -msse2 + runtime
# dispatch). Re-run to refresh after cpp-search updates. See kernel-avx2-reduce-gate.md.
set -e
module load r/4.5.1 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/pjjg18/TreeSearch/lib
BASE=/nobackup/pjjg18/tsmarch; DEPLOY=$BASE/lib; SRC=$BASE/src
mkdir -p $BASE/logs
rm -rf $SRC && git clone --depth 1 -b cpp-search https://github.com/ms609/TreeSearch.git $SRC
echo "cloned cpp-search @ $(cd $SRC && git rev-parse --short HEAD)"
# -ffp-contract=off: keep the AVX2-inline win (the boundary-kill, not FMA) but
# forbid FMA fusion so IW/profile FLOAT scores stay bit-identical to the -msse2
# release -> transparent speedup, safe for IW-tie-sensitive race trajectories.
echo 'PKG_CPPFLAGS = -march=native -ffp-contract=off' > $SRC/src/Makevars
rm -rf $DEPLOY && mkdir -p $DEPLOY
R CMD INSTALL --preclean --no-docs --no-help --library=$DEPLOY $SRC > $BASE/build.log 2>&1 \
  && echo "BUILD OK -> $DEPLOY/TreeSearch" || { echo "BUILD FAILED"; tail -25 $BASE/build.log; exit 1; }
echo "flag on ts_fitch: $(grep 'ts_fitch.cpp' $BASE/build.log | head -1 | grep -oE '\-march=[a-z0-9]+|\-mavx2' | tr '\n' ' ')"
echo "USE: export R_LIBS_USER=$DEPLOY   (or library(TreeSearch, lib.loc='$DEPLOY'))"
echo "### DONE $(date)"
