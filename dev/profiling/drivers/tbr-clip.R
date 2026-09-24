# VTune driver: attribute the per-clip overhead inside the TBR clip loop on the
# FRESH post-fix build.  Pure C++ (one .Call per TBR-to-convergence), so VTune
# sees only ts::tbr_search and its callees.  Zanol2014 = the 2.30x throughput
# outlier (most overhead-bound) -> cleanest attribution.
#
# Goal: rank {compute_insertion_edge_sets (NEW), compute_from_above,
# vroot_cache memcpy, fitch_join_states/fast_hash/dedup (unrooted reroot),
# build_postorder, the bounded scan} by self time.  Recently-added suspects first.
#
# Load the SYMBOLED lib via TS_LIB.  REPS tuned so bare wall is ~3-5 s.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")

fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy <- fc(inapplicable.phyData[["Zanol2014"]])
d <- list(contrast = attr(phy, "contrast"),
          tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
          weight = attr(phy, "weight"), levels = attr(phy, "levels"), nTip = length(phy))

reps <- as.integer(Sys.getenv("REPS", "40"))
t0 <- proc.time()
sc <- integer(reps)
for (i in seq_len(reps)) {
  set.seed(i); start <- RandomTree(phy, root = TRUE)
  edge <- Preorder(RenumberTips(start, names(phy)))[["edge"]]
  set.seed(i)
  res <- TreeSearch:::ts_tbr_diagnostics(edge, d$contrast, d$tip_data, d$weight, d$levels,
           maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, unrooted = TRUE)
  sc[i] <- res$score
}
cat(sprintf("Elapsed: %.1f s  (%d TBR-to-convergence on Zanol2014; score range %d-%d)\n",
            (proc.time() - t0)["elapsed"], reps, min(sc), max(sc)))
