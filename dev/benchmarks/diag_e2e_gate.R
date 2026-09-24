# THE HONEST GATE (advisor): does the plateau fix help END-TO-END, time-matched?
# The oracle (sectorial-from-T0) is null because the 482 basin is across an uphill
# barrier accept-equal can't cross. But the full pipeline has ratchet/drift for
# uphill moves; the question is whether ADDING plateau sector exploration
# (rasStarts>1 + sectorAcceptEqual) helps the full search reach a better score in
# the SAME wall-clock. Time-matched (same maxSeconds) so a win isn't just churn.
#   ON < OFF at matched time => real end-to-end improvement (ship-worthy)
#   ON ~ OFF (or worse)      => plateau is not the closer; gap needs uphill/other
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
dsN   <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
secs  <- as.numeric(Sys.getenv("TS_SECONDS", "120"))
target <- c(Wortley2006 = 482, Zanol2014 = 1262, Zhu2013 = 627, Giles2015 = 671)
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
run <- function(d, seed, plateau) {
  set.seed(seed)
  args <- list(dataset = d, maxSeconds = secs, verbosity = 0L)
  if (plateau) { args$rasStarts <- 3L; args$sectorAcceptEqual <- TRUE }
  r <- suppressWarnings(do.call(MaximizeParsimony, args))
  min(as.double(attr(r, "score")))
}
rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    off <- run(phy, sd, FALSE)
    on  <- run(phy, sd, TRUE)
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd, off = off, on = on,
                                            d = on - off, stringsAsFactors = FALSE)
    cat(sprintf("%-11s s%d | OFF=%.0f ON=%.0f | d=%+.0f | TNT=%s\n",
                nm, sd, off, on, on - off, target[[nm]]))
  }
}
S <- do.call(rbind, rows)
cat("\n== medians (time-matched, full search) ==\n")
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], OFF = median(d$off), ON = median(d$on),
  TNT = target[[d$dataset[1]]])))
print(agg, row.names = FALSE)
