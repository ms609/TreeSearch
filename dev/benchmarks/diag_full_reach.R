# PREMISE RE-CONFIRM (advisor): does TreeSearch's FULL default search ever reach
# TNT's sectorial score? The harness lied about the levers; sanity-check it didn't
# also flatter the headline gap. If full search reaches the target, sectorial-from-T0
# is one weak link others cover (a wall-clock problem). If it plateaus above, the
# sectorial gap is the genuine missing piece.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Zanol2014")), "\\s+")[[1]]
secs <- as.numeric(Sys.getenv("TS_SECONDS", "120"))
target <- c(Wortley2006 = 482, Zanol2014 = 1262, Zhu2013 = 627, Giles2015 = 671)
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  set.seed(1)
  t0 <- proc.time()
  r <- suppressWarnings(MaximizeParsimony(phy, maxSeconds = secs, verbosity = 0L))
  el <- (proc.time() - t0)["elapsed"]
  best <- min(as.double(attr(r, "score")))
  tg <- target[[nm]]
  cat(sprintf("%-11s | full default best=%.0f | TNT target=%s | %s | %.0fs ntrees=%d\n",
              nm, best, ifelse(is.null(tg), "?", tg),
              ifelse(!is.null(tg) && best <= tg, "REACHED", "ABOVE"), el,
              length(r)))
}
