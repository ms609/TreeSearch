# Single-lib timing of an edge-set-heavy EW workload (full auto search,
# Wagner+sector+TBR all call compute_insertion_edge_sets). Run once per lib via
# TS_LIB; a bash loop interleaves libs/rounds and takes medians. Deterministic
# (fixed seed) so any wall delta is the fused-sweep effect, not search-path drift.
lib  <- Sys.getenv("TS_LIB", ".agent-l1")
reps <- as.integer(Sys.getenv("REPS", "10"))
suppressMessages({ library(TreeSearch, lib.loc = lib); library(TreeTools) })
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
for (nm in c("Zhu2013", "Zanol2014")) {
  phy <- fc(inapplicable.phyData[[nm]])
  set.seed(1)
  t0 <- proc.time()
  r <- suppressWarnings(MaximizeParsimony(phy, maxReplicates = reps, nThreads = 1L,
                                          strategy = "auto", verbosity = 0L))
  el <- (proc.time() - t0)["elapsed"]
  cat(sprintf("%s|%s|%.3f|%s\n", lib, nm, el, attr(r, "score")))
}
