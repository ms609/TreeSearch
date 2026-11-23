# Full-search time-matched gate for rasStarts=3 in the AUTO-SELECTED `thorough`
# preset (task #29).  Unlike the rss-only tests, this runs the WHOLE thorough
# pipeline (ratchet/drift/xss/css/rss/fuse interleaved) under a fixed wall-clock
# budget, varying ONLY rasStarts (explicit arg overrides the preset field; all
# other thorough fields unchanged).  This is the decision test before flipping
# thorough's default.  Local wall-clock is INDICATIVE; Hamilton is authoritative.
# Env: TS_LIB, TS_SECONDS, TS_SEEDS, TS_DATASETS.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-sectfix"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
target <- c(Zanol2014 = 1261, Zhu2013 = 624)
dsN   <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Zhu2013")), "\\s+")[[1]]
secs  <- as.integer(Sys.getenv("TS_SECONDS", "60"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2")), "\\s+")[[1]])

run <- function(phy, ras, seed) {
  set.seed(seed)
  r <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
        rasStarts = as.integer(ras), maxSeconds = secs, nThreads = 1L,
        verbosity = 0L))
  min(as.double(attr(r, "score")))
}
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]]); tgt <- target[[nm]]
  cat(sprintf("\n==== %s | target=%d | thorough | budget=%ds ====\n", nm, tgt, secs))
  for (ras in c(1L, 3L)) {
    sc <- vapply(seeds, function(s) run(phy, ras, s), double(1))
    cat(sprintf("  rasStarts=%d -> {%s}  median %+.0f vs target\n",
                ras, paste(sprintf("%.0f", sc), collapse = ","), median(sc) - tgt))
  }
}
