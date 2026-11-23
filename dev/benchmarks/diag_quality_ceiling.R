# QUALITY CEILING: can TreeSearch reach TNT full-xmult's score at all?
#
# headtohead_phase0.csv used strategy="auto" (=thorough for 65-119 tips,
# default for <65) and TS finished 1-4 steps WORSE than TNT xmult:
#   Wortley 483 vs 479 | Zanol 1264 vs 1261 | Zhu 626 vs 624 | Giles 671 vs 670
# But "auto" never tries the strongest preset.  This asks: given the STRONGEST
# TS config (thorough -> intensive -> intensive+intraFuse) and generous time,
# does TS REACH the TNT target, or is it a hard ceiling?
#   REACHED  => gap is preset/tuning/speed (recommend intensive; tune auto).
#   CEILING  => real search-power deficit (need a better component/algorithm).
# nThreads=1 to stay comparable to the headtohead baseline that produced 483/1264.
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
target <- c(Wortley2006 = 479, Eklund2004 = 440, Zanol2014 = 1261,
            Zhu2013 = 624, Giles2015 = 670, Dikow2009 = 1606)

run_arm <- function(phy, seed, secs, extra) {
  set.seed(seed)
  args <- c(list(dataset = phy, maxReplicates = 9999L, maxSeconds = secs,
                 nThreads = 1L, verbosity = 0L), extra)
  r <- tryCatch(suppressWarnings(do.call(MaximizeParsimony, args)),
                error = function(e) { message("  ARM ERROR: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NA_real_)
  min(as.double(attr(r, "score")))
}

arms <- list(
  intensive   = list(strategy = "intensive"),
  plateau     = list(strategy = "intensive", rasStarts = 3L, sectorAcceptEqual = TRUE),
  plateauFuse = list(strategy = "intensive", rasStarts = 3L, sectorAcceptEqual = TRUE,
                     intraFuse = TRUE)
)

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  tgt <- target[[nm]]
  cat(sprintf("\n==== %s (%d tips) | TNT xmult target=%d | %gs x %d seeds, nThreads=1 ====\n",
              nm, NTip(phy), tgt, secs, length(seeds)))
  for (an in names(arms)) {
    sc <- vapply(seeds, function(s) run_arm(phy, s, secs, arms[[an]]), numeric(1))
    best <- suppressWarnings(min(sc, na.rm = TRUE))
    cat(sprintf("  %-14s best=%-5.0f median=%-6.1f all={%s} gap_best=%+.0f%s\n",
                an, best, median(sc, na.rm = TRUE),
                paste(sprintf("%.0f", sc), collapse = ","),
                best - tgt, if (is.finite(best) && best <= tgt) "  <== REACHED" else ""))
  }
}
