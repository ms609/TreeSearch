#!/bin/sh
# Build TreeSearch into the isolated exploration library .agent-softsankoff.
#
# Follows AGENTS.md: tarball into a temp dir, never build in place.  Uses
# CCACHE_DISABLE=1 and --preclean because this branch adds a new header
# (src/ts_soft_sankoff.h); ccache plus a header change is the documented recipe
# for a stale .o compiled against an older struct layout.
set -e

SRC=$(cd "$(dirname "$0")/../.." && pwd)
LIB="$SRC/.agent-softsankoff"
TMPBUILD=$(mktemp -d)

cd "$SRC"
rm -f src/*.o src/*.dll
mkdir -p "$LIB"

cd "$TMPBUILD"
R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC"

CCACHE_DISABLE=1 R CMD INSTALL --preclean --no-docs --no-byte-compile \
  --library="$LIB" "$TMPBUILD"/TreeSearch_*.tar.gz

echo "installed into $LIB"
