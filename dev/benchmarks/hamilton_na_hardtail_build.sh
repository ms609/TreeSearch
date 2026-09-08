#!/bin/bash
# Rebuild TreeSearch from cpp-search for the hard-tail panel.  A fresh build is
# REQUIRED: lib-nacert was built at b9bc14d6, which predates the `effort`
# argument the panel's arms use, so reusing it would fail at the first call.
#SBATCH --job-name=ts-nahard-build
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/nahard_build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/nahard_build_%j.err

module load r/4.5.1
module load gcc/14.2

BRANCH=${TS_BRANCH:-cpp-search}
REPO=/nobackup/$USER/TreeSearch-nahard
DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-nahard
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs
export R_LIBS_USER="$LIB:$DEPLIB"

if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/ms609/TreeSearch.git "$REPO" || { echo "FATAL: clone failed"; exit 1; }
fi
cd "$REPO" || { echo "FATAL: no $REPO"; exit 1; }
git fetch origin "$BRANCH" && (git checkout "$BRANCH" && git reset --hard "origin/$BRANCH")
echo "Git HEAD: $(git log --oneline -1)"

# Clean build: header changes have landed since any earlier build here, and a
# stale .o against an old struct layout segfaults at run time (recorded trap).
rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz

# Fail loudly if `effort` is missing: every arm would silently be the default and
# the panel would report a null for the wrong reason.
Rscript -e '
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  library(TreeSearch)
  if (!"effort" %in% names(formals(MaximizeParsimony))) stop("FATAL: no `effort` -- wrong build")
  data("inapplicable.phyData", package = "TreeSearch")
  d <- inapplicable.phyData[["Vinther2008"]]
  set.seed(1)
  a <- MaximizeParsimony(d, effort = 0L, maxReplicates = 2L, nThreads = 1L, verbosity = 0L)
  Sys.setenv(TS_NA_NOCERTIFY = "1")
  b <- MaximizeParsimony(d, effort = 1L, maxReplicates = 2L, nThreads = 1L, verbosity = 0L)
  if (attr(b, "naDiag")$n_evs_skipped <= 0) stop("FATAL: gate never fired -- wrong build")
  cat("build smoke OK: effort present, gate fires\n")
' || exit 1
exit $rc
