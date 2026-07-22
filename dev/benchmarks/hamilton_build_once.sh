#!/bin/bash
# Build/install TreeSearch ONCE into a shared read-only library, so a job-array
# panel never recompiles. Submit first; chain the array on afterok of this job.
# DISPATCH-UNTESTED template (cell logic validated locally) — smoke in test.q first.
#SBATCH --job-name=ts-build-once
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/build_%j.err

module load r/4.5.1
module load gcc/14.2

LIB=/nobackup/$USER/TreeSearch/lib
REPO=/nobackup/$USER/TreeSearch-a
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs

cd "$REPO" || { echo "FATAL: no $REPO"; exit 1; }
git fetch origin cpp-search && (git checkout cpp-search && git pull --ff-only origin cpp-search \
  || git reset --hard origin/cpp-search)
echo "Git HEAD: $(git log --oneline -1)"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz
echo "INSTALL exit: $rc; version: $(R_LIBS_USER=$LIB Rscript -e 'cat(as.character(packageVersion("TreeSearch")))' 2>/dev/null)"
exit $rc
