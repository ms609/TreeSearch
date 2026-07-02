# T-P5c (Round 7, 2026-07-02): ratchet phase internal wall breakdown.
# Requires a build with the env-gated TS_RATCHET_TIMING chrono instrument in
# ts_ratchet.cpp (measurement-only; NOT in cpp-search — was applied in the
# throwaway worktree TS-ratchet-prof). To re-run: re-apply the instrument
# (three tbr_search calls in ratchet_search + the reject branch), rebuild, point
# TS_LIB at the resulting lib. Splits step1 initial TBR / step2 perturbed-weight
# TBR (scalar scorer, MIXED mode) / step3 restore-to-convergence TBR (flat+x4) /
# reject rebuild. One dataset per process so the cumulative accumulators stay clean.
# RESULT: ratchet = 57% of wall; s1 4% / s2 10% / s3 86% / reject ~0% (Zanol2014,
# thorough, reps3). Perturbed scalar path (s2) is a confirmed non-lever.
Sys.setenv(TS_RATCHET_TIMING = "1")
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }

ds   <- Sys.getenv("TS_DS", "Zanol2014")
reps <- as.integer(Sys.getenv("TS_REPS", "3"))
seed <- as.integer(Sys.getenv("TS_SEED", "1"))
strat <- Sys.getenv("TS_STRAT", "thorough")

phy <- fc(inapplicable.phyData[[ds]])
cat(sprintf("=== %s  nTip=%d  strategy=%s  reps=%d seed=%d ===\n",
            ds, length(phy), strat, reps, seed))
set.seed(seed); t0 <- Sys.time()
r <- suppressWarnings(MaximizeParsimony(
  phy, maxReplicates = reps, nThreads = 1L, strategy = strat, verbosity = 0L))
w <- as.double(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("SCORE=%s  WALL=%.2fs\n", attr(r, "score"), w))
cat("(last [RT] line above = cumulative ratchet sub-phase wall for this run)\n")
