# Before/after for the build_ras_sector directional-edge fix (task #27).
# build_ras_sector (the sector-internal RAS rebuild) fires only when rasStarts>=2,
# so the fix is inert at the default rasStarts=1 (which must therefore be IDENTICAL
# between libs -- a built-in sanity check) and shows only at rasStarts>=2.
#
# rss-only from a shared in-R Wagner start (same on both libs: the build_ras_sector
# change does not touch wagner_tree), seeds aggregated by MaximizeParsimony's own
# multi-start.  Run twice:
#   TS_LIB=.agent-wagsect  Rscript ... > before
#   TS_LIB=.agent-sectfix  Rscript ... > after
# Directional (.agent-sectfix) should be EQUAL-OR-BETTER at rasStarts>=2.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-sectfix"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy <- fitch(inapplicable.phyData[["Zanol2014"]])
target <- 1261

# Shared deterministic Wagner start (identical across libs).
set.seed(7); t0 <- AdditionTree(phy, sequence = sample(seq_along(phy)))
t0len <- TreeLength(t0, phy)
cat(sprintf("lib=%s | Zanol2014 | Wagner T0=%.0f | target=%d\n",
            Sys.getenv("TS_LIB", ".agent-sectfix"), t0len, target))

rss_from <- function(ras) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
        nThreads = 1L, maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L,
        driftCycles = 0L, xssRounds = 0L, cssRounds = 0L, rssRounds = 8L,
        rasStarts = as.integer(ras), wagnerStarts = 1L, fuseInterval = 9999L))
  min(as.double(attr(r, "score")))
}
for (ras in c(1L, 3L, 6L)) {
  sc <- rss_from(ras)
  cat(sprintf("  rss rasStarts=%d -> %.0f  (%+.0f vs T0, %+.0f vs target)\n",
              ras, sc, sc - t0len, sc - target))
}
