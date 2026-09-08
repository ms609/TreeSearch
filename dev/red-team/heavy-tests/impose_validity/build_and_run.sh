#!/usr/bin/env bash
# ---------------------------------------------------------------------
# Build and run the impose_one_pass T-327 guard validity harness.
#
#   bash build_and_run.sh            # n_tip = 4..8 (default)
#   bash build_and_run.sh 4 6        # n_tip = 4..6
#
# Standalone g++ build (no R / Rcpp / Morphy / SIMD): compiles the REAL
# engine src/ts_tree.cpp plus this driver, and #includes the verbatim
# functions-under-test extracted from HEAD.  Invoking g++ directly (not
# `R CMD`) sidesteps this machine's ~/.R/Makevars.win PKG_CXXFLAGS
# clobber — nothing here depends on the R package flags.
# ---------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(git -C "$HERE" rev-parse --show-toplevel)"
SRC="$ROOT/src"

echo "repo root : $ROOT"
echo "extracting verbatim functions-under-test from HEAD ..."
bash "$HERE/extract_funcs.sh" "$ROOT" > "$HERE/extracted_spr.gen.inc"
echo "  -> $(wc -l < "$HERE/extracted_spr.gen.inc") lines generated"

echo "compiling ..."
g++ -std=c++17 -O2 -Wall -Wextra \
    -I"$SRC" -I"$HERE" \
    "$HERE/driver.cpp" "$SRC/ts_tree.cpp" \
    -o "$HERE/driver"

echo "running ..."
echo
"$HERE/driver" "$@"
