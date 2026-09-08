# A/B: xmult-style convergence stop vs the full thorough run.
# Lead lever = consensus_stable_reps (existing param, zero new code): it stops
# when the strict consensus of best-score trees is unchanged for K replicates.
# It is INHERENTLY safe against late score improvements: a new best score
# rebuilds the best-set, changing the consensus hash and resetting the counter,
# so the stability count cannot accumulate while improvements are still arriving.
#
# Verification is on the DELIVERABLE, not just the score:
#   - score must equal the full-run MPT score (no quality loss);
#   - strict-consensus fidelity vs the full run via ClusteringInfoDist (~0);
#   - MPT count retained (early stop must not return a threadbare MPT set).
#
# Also free-simulates a pure score-plateau stop over a K sweep from the full
# run's replicate_scores, to show the (larger, dataset-dependent) K a blunt
# plateau stop would need — contrasting with consensus-stable's safety.
#
# Env: TS_LIB (default .agent-stop), NSEED (default 3).
.libPaths(c(Sys.getenv("TS_LIB", ".agent-stop"), .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools); library(TreeDist) })

nseed <- as.integer(Sys.getenv("NSEED", "3"))
datasets <- c("Wortley2006", "Zanol2014", "Zhu2013", "Giles2015")
target   <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
csK      <- c(3L, 5L)   # consensus_stable_reps arms
data("inapplicable.phyData", package = "TreeSearch")

# Strict consensus as a single phylo (handles the 1-tree case).
.strict <- function(trees) {
  if (inherits(trees, "phylo")) return(trees)
  if (length(trees) == 1L) return(trees[[1]])
  ape::consensus(trees, p = 1)
}
# Smallest score-plateau K that would NOT degrade the final score, simulated
# from the running-min trajectory: walk reps, reset stall on each new min;
# the run stops at the first rep where stall == K. Returns the min K for which
# the achieved (running-min-at-stop) score equals the full final score.
.minSafePlateauK <- function(scores) {
  finalBest <- min(scores)
  runMin <- cummin(scores)
  for (K in seq_len(length(scores))) {
    stall <- 0L; stoppedAt <- length(scores)
    for (i in seq_along(scores)) {
      if (i > 1L && runMin[i] < runMin[i - 1L]) stall <- 0L else stall <- stall + 1L
      if (stall >= K) { stoppedAt <- i; break }
    }
    if (runMin[stoppedAt] <= finalBest) return(K)
  }
  length(scores)
}

runOne <- function(phy, seed, csReps) {
  set.seed(seed)
  t <- system.time(
    r <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
           maxSeconds = 600, nThreads = 1L, verbosity = 0L,
           consensusStableReps = csReps)))
  list(trees = r, score = min(as.double(attr(r, "score"))),
       reps = attr(r, "replicates"),
       lastImp = attr(r, "last_improved_rep"),
       repScores = attr(r, "replicate_scores"),
       consensusStop = isTRUE(attr(r, "consensus_stable")),
       nMPT = length(r), wall = as.double(t["elapsed"]))
}

rows <- list()
for (nm in datasets) {
  m <- PhyDatToMatrix(inapplicable.phyData[[nm]], ambigNA = FALSE); m[m == "-"] <- "?"
  phy <- MatrixToPhyDat(m)
  for (s in seq_len(nseed)) {
    ref <- runOne(phy, s, 0L)                 # full run = reference
    refCons <- .strict(ref$trees)
    safeK <- .minSafePlateauK(ref$repScores)
    for (K in csK) {
      arm <- runOne(phy, s, K)
      armCons <- .strict(arm$trees)
      cid <- tryCatch(
        as.double(ClusteringInfoDist(refCons, armCons, normalize = TRUE)),
        error = function(e) NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = nm, seed = s, csReps = K,
        refScore = ref$score, armScore = arm$score,
        scoreLoss = arm$score - ref$score,
        refReps = ref$reps, armReps = arm$reps,
        refWall = round(ref$wall, 1), armWall = round(arm$wall, 1),
        wallFrac = round(arm$wall / ref$wall, 2),
        refMPT = ref$nMPT, armMPT = arm$nMPT,
        consCID = round(cid, 4), stoppedOnConsensus = arm$consensusStop,
        plateauSafeK = safeK)
      cat(sprintf(paste0("%-12s s%d cs%d: score %.0f vs %.0f (loss %+.0f) | ",
                         "reps %d->%d  wall %.1f->%.1fs (%.0f%%) | ",
                         "MPT %d->%d  consCID %.3f | plateauSafeK=%d %s\n"),
                  nm, s, K, ref$score, arm$score, arm$score - ref$score,
                  ref$reps, arm$reps, ref$wall, arm$wall, 100 * arm$wall / ref$wall,
                  ref$nMPT, arm$nMPT, ifelse(is.na(cid), -1, cid), safeK,
                  ifelse(arm$consensusStop, "[cons-stop]", "[other-stop]")))
    }
  }
}
df <- do.call(rbind, rows)
write.csv(df, file.path(Sys.getenv("OUTDIR", "dev/benchmarks"),
                        "convergence_ab.csv"), row.names = FALSE)
cat("\n=== median by dataset x csReps ===\n")
agg <- aggregate(cbind(scoreLoss, armReps, wallFrac, armMPT, consCID, plateauSafeK) ~
                   dataset + csReps, df, median)
print(agg[order(agg$dataset, agg$csReps), ], row.names = FALSE)
cat("\nGo if: scoreLoss==0 everywhere, consCID~0, armMPT not collapsed, wallFrac<<1.\n")
cat("plateauSafeK shows the (larger) K a blunt score-plateau would need.\n")
