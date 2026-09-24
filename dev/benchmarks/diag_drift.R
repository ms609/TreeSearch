# DRIFT LEG: TNT's sectsch plateaus at the TS ceiling on Wortley (480); the final
# 480->479 is xmult's DRIFT (+fuse), which every TreeSearch preset disables
# (driftCycles=0, audit D3).  Does enabling TNT-faithful drift let TS reach the
# xmult target?  Wortley 479 (sectsch-null, drift-only) is the clean case; Zanol
# 1261 is mixed (sectorial + drift).  nThreads=1, best-of-N over seeds.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
secs   <- as.numeric(Sys.getenv("TS_SECONDS", "60"))
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006")), "\\s+")[[1]]
seeds  <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5")), "\\s+")[[1]])
target <- c(Wortley2006 = 479, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)

run <- function(phy, seed, extra) {
  set.seed(seed)
  args <- c(list(dataset = phy, maxReplicates = 9999L, maxSeconds = secs,
                 nThreads = 1L, verbosity = 0L), extra)
  r <- tryCatch(suppressWarnings(do.call(MaximizeParsimony, args)),
                error = function(e) { message("ERR: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NA_real_)
  min(as.double(attr(r, "score")))
}

arms <- list(
  intensive          = list(strategy = "intensive"),                                  # drift=0 baseline
  drift30            = list(strategy = "intensive", driftCycles = 30L),
  drift30_fuse       = list(strategy = "intensive", driftCycles = 30L, intraFuse = TRUE)
)
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  tgt <- target[[nm]]
  cat(sprintf("\n==== %s (%d tips) | xmult target=%d | %gs x %d seeds ====\n",
              nm, NTip(phy), tgt, secs, length(seeds)))
  for (an in names(arms)) {
    sc <- vapply(seeds, function(s) run(phy, s, arms[[an]]), numeric(1))
    best <- suppressWarnings(min(sc, na.rm = TRUE))
    cat(sprintf("  %-14s best=%-5.0f median=%-6.1f all={%s} gap=%+.0f%s\n",
                an, best, median(sc, na.rm = TRUE), paste(sprintf("%.0f", sc), collapse = ","),
                best - tgt, if (is.finite(best) && best <= tgt) "  <== REACHED" else ""))
  }
}
