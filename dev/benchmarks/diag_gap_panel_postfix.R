# Gap-panel re-measurement AFTER the Wagner + TBR-vroot + build_ras_sector
# directional-cost fixes.  The 2026-06-16 plan concluded the EW-Fitch score gap
# was at a "landscape/escape-bound floor" (+1.5..+3.5 on the hard datasets) with
# a "competitive per-candidate" kernel -- but that predated finding the
# union-of-finals cost bug.  This re-runs the full `thorough` pipeline at a fixed
# budget on the hard panel to test whether the bug fix closed the score gap.
# Env: TS_LIB (default .agent-sectfix), TS_SECONDS, TS_SEEDS, TS_DATASETS.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-sectfix"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
# Known EW-Fitch MPT / TNT targets (apples-to-apples, -> ?).
target <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
dsN   <- strsplit(trimws(Sys.getenv("TS_DATASETS",
           "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
secs  <- as.integer(Sys.getenv("TS_SECONDS", "60"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])

cat(sprintf("Gap panel | thorough | %ds | seeds {%s} | lib %s\n",
            secs, paste(seeds, collapse=","), Sys.getenv("TS_LIB", ".agent-sectfix")))
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]]); tgt <- target[[nm]]
  sc <- vapply(seeds, function(s) {
    set.seed(s)
    r <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
           maxSeconds = secs, nThreads = 1L, verbosity = 0L))
    min(as.double(attr(r, "score")))
  }, double(1))
  cat(sprintf("  %-12s target=%4d  scores {%s}  min %+.0f  median %+.0f\n",
              nm, tgt, paste(sprintf("%.0f", sc), collapse=","),
              min(sc) - tgt, median(sc) - tgt))
}
