#!/bin/bash
#SBATCH --job-name=2island-build
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/2island_build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/2island_build_%j.err
#
# Build step for the two-island "thorough" tuning sweep. Builds TreeSearch from
# origin/cpp-search into a shared lib, then the array job (hamilton_two_island_
# sweep.sh) runs against it. Submit with a dependency:
#   bid=$(sbatch --parsable hamilton_two_island_build.sh)
#   sbatch --dependency=afterok:$bid hamilton_two_island_sweep.sh
#
# PREREQUISITE: builds origin/cpp-search-2island (= cpp-search + the two
# benchmark commits 356fc8d0 + b9164a77, which carry the vendored fixture and
# this driver). Override the ref with TS_BRANCH if you graft them elsewhere.

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-2island
LIB=/nobackup/$USER/TreeSearch/lib2island
# Branch carrying the vendored two-island fixture + sweep driver. Defaults to
# cpp-search-2island (cpp-search + the two benchmark commits, grafted cleanly).
BRANCH=${TS_BRANCH:-cpp-search-2island}
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs
export R_LIBS="$LIB:${R_LIBS}"

echo "=== two-island sweep build === $(date) | node $(hostname)"

# Dependencies (login node normally has these; install any missing).
Rscript --no-save -e "
  lib <- '$LIB'; .libPaths(c(lib, .libPaths()))
  pkgs <- c('abind','ape','cli','colorspace','fastmatch','Rdpack','RcppParallel','TreeDist','TreeTools','PlotTools')
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) install.packages(need, lib = lib,
        repos = 'https://cloud.r-project.org', dependencies = NA, quiet = TRUE)
"

if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/ms609/TreeSearch.git "$REPO"
fi
cd "$REPO" || exit 1
git fetch origin "$BRANCH" && git reset --hard "origin/$BRANCH"
echo "Git branch: $BRANCH | HEAD: $(git log --oneline -1)"

# Confirm the vendored fixture made it (else the sweep's Part A cannot run).
for f in dev/benchmarks/data/zhu2013_2island.nex \
         dev/benchmarks/data/zhu2013_2island_island2.tre \
         dev/benchmarks/data/zhu2013_2island_main.tre; do
  [ -f "$f" ] || { echo "FATAL: missing fixture $f -- push cpp-search first"; exit 2; }
done

echo "Building TreeSearch..."
rm -f src/*.o src/*.so src/*.dll
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz; rc=$?
rm -f TreeSearch_*.tar.gz
if [ $rc -ne 0 ]; then echo "FATAL: install failed"; exit 1; fi

touch "$LIB/.build_done"
echo "=== build complete === $(date)"
