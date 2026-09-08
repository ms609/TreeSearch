# Scoring-approximation Δ-probe driver (T-F1 / task #53).
# Forces prune_reinsert (expand_and_reinsert) on a real heavy-multistate dataset
# and lets the -DTS_SCOREAPPROX_PROBE instrumentation tally, per tip placement:
#   Delta = exact_cost(E_bounded) - min_E exact_cost(E)
# i.e. how much GREEDY insertion cost the undercounting union-of-finals scorer
# loses vs the exact directional edge_set scorer.  Cumulative [SCOREAPPROX] line
# is printed to stderr after each MaximizeParsimony() call; the LAST line is the
# grand total for this process.  Observational only: production still inserts at
# the _bounded choice, so the search trajectory is unchanged.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy  <- fc(inapplicable.phyData[[Sys.getenv("TS_DS", "Zanol2014")]])
cyc  <- as.integer(Sys.getenv("TS_PRC",  "30"))
drop <- as.double (Sys.getenv("TS_DROP", "0.10"))
sel  <- as.integer(Sys.getenv("TS_SEL",  "0"))
reps <- as.integer(Sys.getenv("TS_REPS", "4"))
cat(sprintf("DS=%s tips=%d | pruneReinsertCycles=%d drop=%.2f sel=%d seeds=%d\n",
            Sys.getenv("TS_DS", "Zanol2014"), NTip(phy), cyc, drop, sel, reps))
t0 <- proc.time()
for (i in seq_len(reps)) {
  set.seed(1000L + i)
  invisible(suppressWarnings(MaximizeParsimony(
    phy, maxReplicates = 1L, nThreads = 1L, strategy = "default", verbosity = 0L,
    ratchetCycles = 0L,
    pruneReinsertCycles = cyc, pruneReinsertDrop = drop, pruneReinsertSelection = sel)))
}
cat("Elapsed:", round((proc.time() - t0)["elapsed"], 1), "s\n")
