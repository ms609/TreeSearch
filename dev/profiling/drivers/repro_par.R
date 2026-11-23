# Tight reproducer for the nThreads=2 Vinther2008 (NA path) crash.
# Runs the exact failing scenario N times; prints before each so an abnormal
# exit pinpoints the iteration. Lib via TS_LIB. Reports survivors + scores.
lib <- Sys.getenv("TS_LIB", ".agent-l1")
n   <- as.integer(Sys.getenv("N", "6"))
suppressMessages({ library(TreeSearch, lib.loc = lib); library(TreeTools) })
data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]
cat(sprintf(">>> lib=%s  n=%d  (Vinther2008, NA, nThreads=2)\n", lib, n)); flush.console()
ok <- 0L
for (i in seq_len(n)) {
  cat(sprintf("  iter %d ... ", i)); flush.console()
  set.seed(6274 + i)
  r <- suppressWarnings(MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                                          nThreads = 2L, verbosity = 0L))
  cat(sprintf("score=%s class=%s\n", attr(r, "score"), class(r)[1])); flush.console()
  ok <- ok + 1L
}
cat(sprintf(">>> SURVIVED %d/%d\n", ok, n))
