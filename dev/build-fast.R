#!/usr/bin/env Rscript
# Fast dev build for C++-only iteration.
#
# Incremental -O2 compile (only changed translation units, via ccache + parallel
# make from ~/.R/Makevars.win), then HOT-SWAP the freshly built DLL into a target
# install library so benchmarks/tests pick it up WITHOUT a full R CMD INSTALL.
#
# Use this for C++-only edits. For changes to R/, roxygen, or [[Rcpp::export]]
# SIGNATURES, do a full `R CMD INSTALL` (compileAttributes is run here to catch
# export changes, but a signature change still needs the R wrapper reinstalled).
#
# debug = FALSE => -O2, so timing / candidate-throughput / profiling stay valid.
# (Use compile_dll(debug=TRUE) -> -O0 ONLY for logic/correctness loops, never timing.)
#
# Usage: Rscript dev/build-fast.R [target_lib=.agent-p0]
#   Then run benchmarks against that lib (lib.loc / TS_LIB = target_lib).
#   Do NOT have an R session holding the target DLL open (Windows file lock).

args <- commandArgs(trailingOnly = TRUE)
lib  <- if (length(args) >= 1L) args[[1]] else ".agent-p0"

t0 <- Sys.time()
Rcpp::compileAttributes(".")              # guards the stale-RcppExports trap; no-op if unchanged
pkgbuild::compile_dll(".", debug = FALSE) # incremental, -O2; recompiles only changed TUs
build_s <- as.double(difftime(Sys.time(), t0, units = "secs"))

dll <- Sys.glob(file.path("src", "*.dll"))
dst <- file.path(lib, "TreeSearch", "libs", "x64", "TreeSearch.dll")
if (length(dll) == 1L && dir.exists(dirname(dst))) {
  ok <- file.copy(dll, dst, overwrite = TRUE)
  cat(sprintf("Hot-swapped %s -> %s  [%s]\n", dll, dst,
              if (ok) "ok" else "FAILED (DLL locked? close R sessions using this lib)"))
} else {
  cat(sprintf("No hot-swap: dll matches=%d, target dir exists=%s. Run a full install into '%s' first.\n",
              length(dll), dir.exists(dirname(dst)), lib))
}
cat(sprintf("build-fast: %.1fs\n", build_s))
