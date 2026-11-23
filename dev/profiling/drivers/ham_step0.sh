#!/bin/bash
#SBATCH --job-name=ts-step0-64bit
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=1:30:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/step0_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/step0_%j.err
set -e
module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LIB=/nobackup/$USER/TreeSearch/lib
REPO=/nobackup/$USER/TreeSearch-t29
export R_LIBS_USER=$LIB
mkdir -p /nobackup/$USER/TreeSearch/logs

cd "$REPO"
# Repo already fast-forwarded to origin/cpp-search@78b74147 on the login node
# (unrooted-default 25e35be7 present); compute nodes may lack outbound network.
echo "BUILD HEAD: $(git log --oneline -1)"
git merge-base --is-ancestor 25e35be7 HEAD && echo "unrooted-default: YES" || echo "unrooted-default: NO"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
echo "INSTALL exit: $?"
rm -f TreeSearch_*.tar.gz

# 64-bit gold-standard like-vs-like (framing.R = fitch-only rate/gapB decomposition)
export TS_LIB=$LIB
export TNT_EXE=/nobackup/$USER/TreeSearch/tnt/TNT-bin/tnt
export TS_DATASETS="Wortley2006 Wills2012 Zanol2014 Zhu2013 Giles2015"
export TS_SEEDS="1 2 3"
export TS_SECONDS=300
export OUT_CSV=/nobackup/$USER/TreeSearch/framing_64bit.csv
Rscript /nobackup/$USER/TreeSearch/framing.R
echo "DONE step0"
