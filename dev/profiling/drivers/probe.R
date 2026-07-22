# Diagnostic probe: time each MaximizeParsimony run against ONE library,
# printing per-run progress so a hang can be localized.  NOT the gate.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Zhu2013 Zanol2014")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1")), "\\s+")[[1]])
cat(sprintf("LIB=%s\n", Sys.getenv("TS_LIB")))
for (nm in dsN) {
  raw <- inapplicable.phyData[[nm]]; fit <- fc(raw)
  for (obj in c("fitch", "raw")) {
    phy <- if (obj == "fitch") fit else raw
    for (sd in seeds) {
      cat(sprintf("[start] %s/%s/seed%d ... ", nm, obj, sd)); flush.console()
      set.seed(sd); t0 <- Sys.time()
      r <- suppressWarnings(MaximizeParsimony(
        phy, maxReplicates = 3L, nThreads = 1L, strategy = "auto", verbosity = 0L))
      w <- as.double(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("done score=%s cand=%s wall=%.2fs\n",
                  attr(r, "score"), attr(r, "candidates_evaluated"), w)); flush.console()
    }
  }
}
cat("PROBE COMPLETE\n")
