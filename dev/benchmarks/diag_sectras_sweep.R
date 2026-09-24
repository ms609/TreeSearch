# How far does raising rasStarts close the sectorial gap, and at what cost?
# rss-only from a shared in-R Wagner T0, per dataset, rasStarts in {1,3,6}, with
# wall-clock. Scores are bitness-independent (vs the hardcoded TNT targets);
# relative timing shows the rasStarts cost multiplier. Informs whether a
# TNT-faithful rasStarts (TNT uses 3) should be the sectorial default/preset.
# Env: TS_LIB (default .agent-sectfix), TS_DATASETS.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-sectfix"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
target <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Zhu2013 Wortley2006")), "\\s+")[[1]]

rss_from <- function(phy, t0, ras) {
  set.seed(1)
  t <- system.time(r <- suppressWarnings(MaximizeParsimony(phy, tree = t0,
        maxReplicates = 1L, nThreads = 1L, maxSeconds = 0, verbosity = 0L,
        ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
        rssRounds = 8L, rasStarts = as.integer(ras), wagnerStarts = 1L,
        fuseInterval = 9999L)))
  list(score = min(as.double(attr(r, "score"))), secs = as.double(t["elapsed"]))
}

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  set.seed(7); t0 <- AdditionTree(phy, sequence = sample(seq_along(phy)))
  t0len <- TreeLength(t0, phy); tgt <- target[[nm]]
  cat(sprintf("\n==== %s | Wagner T0=%.0f | target=%d ====\n", nm, t0len, tgt))
  for (ras in c(1L, 3L, 6L)) {
    o <- rss_from(phy, t0, ras)
    cat(sprintf("  rasStarts=%d -> %.0f  (%+.0f vs target)  [%.1fs]\n",
                ras, o$score, o$score - tgt, o$secs))
  }
}
