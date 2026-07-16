#!/bin/bash
# Build TreeSearch @ claude/dss-probe ONCE into a dedicated lib for the DSS size-axis
# probe. Uses a DEDICATED repo dir + branch (not the shared TreeSearch-a/cpp-search
# checkout) so it never collides with concurrent Hamilton jobs. Submit first; chain
# the array on afterok of this job.
#SBATCH --job-name=ts-dss-build
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/dssbuild_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/dssbuild_%j.err

module load r/4.5.1
module load gcc/14.2

LIB=/nobackup/$USER/TreeSearch/dsslib
REPO=/nobackup/$USER/TreeSearch-dss
SRC=/nobackup/$USER/TreeSearch-a           # existing clone, source of origin URL
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs

# Dedicated checkout of the probe branch (clone once from origin, then fetch).
if [ ! -d "$REPO/.git" ]; then
  ORIGIN=$(git -C "$SRC" remote get-url origin)
  echo "Cloning $ORIGIN -> $REPO"
  git clone "$ORIGIN" "$REPO" || { echo "FATAL: clone failed"; exit 1; }
fi
cd "$REPO" || { echo "FATAL: no $REPO"; exit 1; }
git fetch origin claude/dss-probe && git checkout -B claude/dss-probe origin/claude/dss-probe \
  && git reset --hard origin/claude/dss-probe
echo "Git HEAD: $(git log --oneline -1)"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz
echo "INSTALL exit: $rc; version: $(R_LIBS_USER=$LIB Rscript -e 'cat(as.character(packageVersion("TreeSearch")))' 2>/dev/null)"
exit $rc
