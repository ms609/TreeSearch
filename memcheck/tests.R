# TARGETED ASan driver (diagnostic branch claude/asan-na-incr-crash).
#
# Pins the Windows-only crash in test-ts-na-incremental.R Test 11: the NA
# incremental exact_verify rescore (default-on since #254) crashes after the
# audit prints "5658/5658 byte-matched" — abnormal termination, not a testthat
# failure (truncated .Rout.fail). Local Windows could not reproduce it (heap-
# layout-dependent); gcc-ASan instruments every access and will report the
# offending file:line if it is a generic OOB / use-after-free.
#
# Runs BOTH the shipping path (production incremental) and the audit path on the
# exact crash config (Vinther2008, IW concavity=10, as.phylo(1)). Which one ASan
# flags decides whether the fix is in the shared kernel or audit-only scaffolding.
suppressMessages({
  library("TreeSearch")
  library("TreeTools")
})
data("inapplicable.phyData", package = "TreeSearch")

mk <- function(phy) list(
  contrast = attr(phy, "contrast"),
  tip_data = matrix(unlist(phy, use.names = FALSE),
                    nrow = length(phy), byrow = TRUE),
  weight   = attr(phy, "weight"),
  levels   = attr(phy, "levels"))

dataset  <- inapplicable.phyData[["Vinther2008"]]
ds       <- mk(dataset)
n_tip    <- length(dataset)
minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))
start    <- as.phylo(1, n_tip)

climb <- function() {
  TreeSearch:::ts_tbr_search(start$edge, ds$contrast, ds$tip_data, ds$weight,
                             ds$levels, maxHits = 1L, min_steps = minSteps,
                             concavity = 10)
}

cat(">>> PRODUCTION incremental (shipping path) <<<\n"); flush.console()
Sys.unsetenv("TS_NA_INCR_AUDIT"); Sys.unsetenv("TS_NA_NOINCR")
for (i in 1:3) {
  set.seed(123); r <- climb()
  cat(sprintf("  production iter %d OK score=%.6g\n", i, r$score)); flush.console()
}

cat(">>> AUDIT path (reproduces CI Test 11) <<<\n"); flush.console()
Sys.setenv(TS_NA_INCR_AUDIT = "1")
for (i in 1:2) {
  set.seed(123); r <- climb()
  cat(sprintf("  audit iter %d OK score=%.6g\n", i, r$score)); flush.console()
}
Sys.unsetenv("TS_NA_INCR_AUDIT")

cat("DRIVER COMPLETE (no ASan abort)\n")
