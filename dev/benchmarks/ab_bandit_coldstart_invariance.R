#!/usr/bin/env Rscript
# Invariance check for the warm-start bandit fix: with no `tree = ` supplied,
# every replicate builds its own start, `user_started` is always false, and the
# guard is unchanged — so `before` and `after` must agree EXACTLY on both score
# and the per-arm attempt vector, for every seed.
#
# This is the empirical backing for the NEWS claim that searches using
# `adaptiveStart` alone are unaffected.  Any mismatch falsifies it.
#
# Usage: Rscript dev/benchmarks/ab_bandit_coldstart_invariance.R

beforeLib <- ".ab-bandit-before"
afterLib  <- ".ab-bandit-after"
datasets  <- c("Vinther2008", "Wortley2006")
seeds     <- 1:3
nRep      <- 10L

runCold <- function(lib, dsName, seed) {
  tmp <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', lib),
    'suppressMessages(library(TreeSearch))',
    sprintf('ds <- inapplicable.phyData[["%s"]]', dsName),
    sprintf('set.seed(%d)', seed),
    sprintf('r <- MaximizeParsimony(ds, maxReplicates = %dL, targetHits = 99L,',
            nRep),
    '  adaptiveStart = TRUE, verbosity = 0L, nThreads = 1L)',
    'cat(sprintf("%g|%s\\n", attr(r, "score"),',
    '  paste(attr(r, "strategy_diagnostics")$attempts, collapse = ",")))'
  ), tmp)
  # stderr passes through to the console: a crashed child must not look like a
  # quiet mismatch.
  out <- system2("Rscript", c("--no-save", tmp), stdout = TRUE, stderr = "")
  unlink(tmp)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    return(sprintf("<child exited %d>", status))
  }
  if (length(out) == 0L) return("<no output>")
  trimws(tail(out, 1))
}

nMismatch <- 0L
for (dsName in datasets) {
  for (seed in seeds) {
    b <- runCold(beforeLib, dsName, seed)
    a <- runCold(afterLib,  dsName, seed)
    same <- identical(b, a)
    if (!same) nMismatch <- nMismatch + 1L
    cat(sprintf("%-14s seed %d  %s\n    before %s\n    after  %s\n",
                dsName, seed, if (same) "IDENTICAL" else "*** DIFFERS ***",
                b, a))
  }
}

cat(sprintf("\n%d/%d cells identical; %d mismatch(es)\n",
            length(datasets) * length(seeds) - nMismatch,
            length(datasets) * length(seeds), nMismatch))
if (nMismatch > 0L) {
  stop("Cold-start invariance VIOLATED — the fix is not confined to warm starts")
}
cat("Cold-start invariance holds.\n")
