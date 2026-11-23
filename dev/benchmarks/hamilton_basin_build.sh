#!/bin/bash
#SBATCH --job-name=ts_basin_build
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/basin_build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/basin_build_%j.err
#
# Build TreeSearch ONCE into an isolated per-measurement library, so the basin
# array never recompiles. Deps must already be installed into $LIB on the login
# node (see the dispatch notes). Chain the array on afterok of this job.
set -euo pipefail

module load r/4.5.1 gcc/14.2

: "${REPO:=/nobackup/$USER/TreeSearch-basin}"
: "${LIB:=/nobackup/$USER/TreeSearch-basin-lib}"
# Dependency base: a lib already holding ape/Rcpp/TreeTools/TreeDist/... so we
# compile only TreeSearch itself into the isolated $LIB.
: "${DEPS_LIB:=/nobackup/$USER/TreeSearch/lib}"
export R_LIBS_USER="$LIB:$DEPS_LIB"
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs

cd "$REPO"
echo "Git HEAD: $(git log --oneline -1)"
rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz
echo "INSTALL exit: $rc; version: $(Rscript -e 'cat(as.character(packageVersion("TreeSearch")))' 2>/dev/null)"
exit $rc
