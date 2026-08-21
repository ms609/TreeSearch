#!/bin/bash
# T-364/T-370 three-arm constrained-search battery: build ONE arm into its own
# library.  Submit three times with --export=ARM=1|2|3.
#
#   ARM 1 = bbcca1ba  pre-fix (parent of 7685bf07)
#   ARM 2 = 7685bf07  complement enforcement only, NO reroot  <- T-384 exposed
#   ARM 3 = 796a29d3  complement + reroot at tip 0             <- merged state
#
# Each arm gets its OWN clone, so every build starts with zero object files:
# no stale .o can survive a header-layout change (memory stale-object-abi-gotcha).
#SBATCH --job-name=t364-build
#SBATCH -p shared
#SBATCH -n 8
#SBATCH --mem=12G
#SBATCH --time=1:30:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t364build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t364build_%j.err

module load r/4.5.1

ARM=${ARM:?set ARM=1|2|3}
case "$ARM" in
  1) SHA=bbcca1ba ;;
  2) SHA=7685bf07 ;;
  3) SHA=796a29d3 ;;
  *) echo "FATAL: bad ARM=$ARM"; exit 2 ;;
esac

BASE=/nobackup/$USER/TreeSearch
REPO=$BASE-t364arm$ARM
LIB=$BASE/t364-lib$ARM
DEPLIB=$BASE/lib
mkdir -p "$LIB" "$BASE/logs"

if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/ms609/TreeSearch.git "$REPO" || exit 1
fi
cd "$REPO" || exit 1
git fetch origin cpp-search || exit 1
git checkout --detach "$SHA" || exit 1
echo "ARM $ARM HEAD: $(git log --oneline -1)"
echo "ARM $ARM expected SHA prefix: $SHA ; actual: $(git rev-parse --short=8 HEAD)"

rm -f src/*.o src/*.so
export MAKEFLAGS="-j8"
R_LIBS="$DEPLIB" R CMD INSTALL --no-docs --library="$LIB" .
rc=$?
echo "ARM $ARM INSTALL exit: $rc"
R_LIBS="$LIB:$DEPLIB" Rscript -e 'cat("version:", as.character(packageVersion("TreeSearch")), "\n")'
exit $rc
