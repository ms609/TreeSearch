#!/bin/bash
#SBATCH --job-name=ts-step0-tntfix
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/step0rr_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/step0rr_%j.err
set -e
module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# THE FIX (verified 2026-06-19): the earlier "curses needs TERM" theory was WRONG.
# 64-bit TNT starts fine headless; it died because `tnt run.run;` makes TNT try to
# EXECUTE the filename as a command and `proc` never runs.  framing.R now feeds the
# script via STDIN (the canonical non-interactive mode); no curses is initialised
# when stdin is not a TTY, so TERM=dumb suffices.  Results land in the `log` file.
export TERM=dumb
LIB=/nobackup/$USER/TreeSearch/lib
export R_LIBS_USER=$LIB
# REUSE the lib installed by job 17528864 (TS rates already validated); NO rebuild.
Rscript -e 'cat("lib check: TreeSearch", as.character(packageVersion("TreeSearch")), "\n")'

export TS_LIB=$LIB
export TNT_EXE=/nobackup/$USER/TreeSearch/tnt/TNT-bin/tnt
export TS_DATASETS="Wortley2006 Wills2012 Zanol2014 Zhu2013 Giles2015"
export TS_SEEDS="1 2 3"
export TS_SECONDS=300
export OUT_CSV=/nobackup/$USER/TreeSearch/framing_64bit_log.csv
Rscript /nobackup/$USER/TreeSearch/framing.R
echo "DONE step0-rerun"
