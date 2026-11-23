# TIME-MATCHED rss: rasStarts=1 (many shallow rounds) vs 3 (fewer deeper rounds)
# under an IDENTICAL wall-clock budget. The unbounded sweep (diag_sectras_sweep.R)
# showed ras=3 reaches +1 vs ras=1's +7/+8 when rss runs to completion -- but ras=3
# costs ~3-5x/sector, so the real question for a preset change is whether it still
# wins when TIME is the constraint. rssRounds set high so maxSeconds is the bound.
# Local wall-clock is only indicative (Hamilton is authoritative); the score
# comparison at matched budget is the signal. Env: TS_LIB, TS_SECONDS, TS_SEEDS.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-sectfix"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
target <- c(Zanol2014 = 1261, Zhu2013 = 624)
dsN  <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Zhu2013")), "\\s+")[[1]]
secs <- as.integer(Sys.getenv("TS_SECONDS", "30"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2")), "\\s+")[[1]])

rss_timed <- function(phy, t0, ras, seed) {
  set.seed(seed)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
        nThreads = 1L, maxSeconds = secs, verbosity = 0L, ratchetCycles = 0L,
        driftCycles = 0L, xssRounds = 0L, cssRounds = 0L, rssRounds = 50L,
        rasStarts = as.integer(ras), wagnerStarts = 1L, fuseInterval = 9999L))
  min(as.double(attr(r, "score")))
}
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  set.seed(7); t0 <- AdditionTree(phy, sequence = sample(seq_along(phy)))
  tgt <- target[[nm]]
  cat(sprintf("\n==== %s | target=%d | budget=%ds ====\n", nm, tgt, secs))
  for (ras in c(1L, 3L)) {
    sc <- vapply(seeds, function(s) rss_timed(phy, t0, ras, s), double(1))
    cat(sprintf("  rasStarts=%d -> scores {%s}  median %+.0f vs target\n",
                ras, paste(sprintf("%.0f", sc), collapse = ","),
                median(sc) - tgt))
  }
}
