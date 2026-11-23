# Convergence-tail diagnostic: is there a recoverable wall-clock tail between
# time-to-first-MPT and the stop replicate?  TS self-terminates on targetHits
# (~nTip/5 rediscoveries); TNT xmult stops the moment the score stops improving.
# Because best_score is monotonic, `last_improved_rep` IS time-to-first-MPT.
# The gap (replicates - last_improved_rep) is the tail an xmult-style stop cuts.
#
# Also records the MAX gap between consecutive score improvements on the path to
# the optimum (from replicate_scores), which sets the floor for any plateau-K.
#
# Env: TS_LIB (default .agent-stop), NSEED (default 3).
.libPaths(c(Sys.getenv("TS_LIB", ".agent-stop"), .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools) })

nseed <- as.integer(Sys.getenv("NSEED", "3"))
datasets <- c("Wortley2006", "Zanol2014", "Zhu2013", "Giles2015")
target   <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)

data("inapplicable.phyData", package = "TreeSearch")

# Running-minimum of the per-replicate best score, and the replicate at which
# each new minimum (improvement) first appeared; returns the max gap between
# consecutive improvements (patience a plateau-stop must survive).
.maxImproveGap <- function(scores) {
  if (length(scores) < 2L) return(0L)
  runMin <- cummin(scores)
  improveReps <- which(c(TRUE, diff(runMin) < 0))   # reps that lowered the min
  if (length(improveReps) < 2L) return(0L)
  max(diff(improveReps))
}

rows <- list()
for (nm in datasets) {
  m <- PhyDatToMatrix(inapplicable.phyData[[nm]], ambigNA = FALSE)
  m[m == "-"] <- "?"
  phy <- MatrixToPhyDat(m)
  nTip <- length(phy)
  for (s in seq_len(nseed)) {
    set.seed(s)
    t <- system.time(
      r <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
             maxSeconds = 600, nThreads = 1L, verbosity = 0L)))
    reps      <- attr(r, "replicates")
    lastImp   <- attr(r, "last_improved_rep")
    repScores <- attr(r, "replicate_scores")
    stopReason <- if (isTRUE(attr(r, "consensus_stable"))) "consensus"
                  else if (isTRUE(attr(r, "perturb_stop"))) "perturb"
                  else "targetHits/max"
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = nm, nTip = nTip, target = target[[nm]], seed = s,
      score = min(as.double(attr(r, "score"))),
      reps = reps, lastImproveRep = lastImp,
      tailReps = reps - lastImp,                 # recoverable replicates
      maxImproveGap = .maxImproveGap(repScores), # floor for plateau-K
      stopReason = stopReason,
      wall_s = round(as.double(t["elapsed"]), 1))
    cat(sprintf("%-12s s%d: score %.0f (%+.0f) | reps %d, last-improve %d, tail %d, maxGap %d | %s | %.1fs\n",
                nm, s, min(as.double(attr(r, "score"))),
                min(as.double(attr(r, "score"))) - target[[nm]],
                reps, lastImp, reps - lastImp, .maxImproveGap(repScores),
                stopReason, t["elapsed"]))
  }
}
df <- do.call(rbind, rows)
outdir <- Sys.getenv("OUTDIR", "dev/benchmarks")
write.csv(df, file.path(outdir, "convergence_tail.csv"), row.names = FALSE)
cat("\n=== median by dataset ===\n")
agg <- aggregate(cbind(reps, lastImproveRep, tailReps, maxImproveGap, wall_s) ~ dataset,
                 df, median)
print(agg, row.names = FALSE)
cat("\nIf tailReps >> 0 and maxImproveGap is small, an xmult-style plateau stop ",
    "recovers (tailReps/reps) of the wall.\n", sep = "")
