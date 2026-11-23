# VTune driver: FULL EW MaximizeParsimony (the real-world / mission workload),
# to localize the per-ITERATION deficit vs TNT. ~90% of full-EW time is TBR
# candidate evaluation (ratchet+sectorial+TBR), so flat hotspots here attribute
# the per-candidate cost: the per-clip RECOMPUTE cluster (compute_insertion_edge_sets
# + compute_from_above + vroot_cache rebuild + build_postorder) vs irreducible
# core Fitch scoring (fitch_indirect_length_bounded, any_hit_reduce).
# Loads the symboled lib named by TS_LIB. Target ~5s bare (VTune adds 5-20x).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy   <- fc(inapplicable.phyData[[Sys.getenv("TS_DS", "Zanol2014")]])
reps  <- as.integer(Sys.getenv("TS_REPS", "3"))
iters <- as.integer(Sys.getenv("TS_ITERS", "2"))
t0 <- proc.time()
for (i in seq_len(iters)) {
  set.seed(i)
  invisible(suppressWarnings(MaximizeParsimony(
    phy, maxReplicates = reps, nThreads = 1L, strategy = "auto", verbosity = 0L)))
}
cat("Elapsed:", round((proc.time() - t0)["elapsed"], 2), "s\n")
