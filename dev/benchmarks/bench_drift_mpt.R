#!/usr/bin/env Rscript
# T-254: Drift MPT diversity experiment
#
# Compare pool size, MPT count, and topological diversity between
# driftCycles=0 and driftCycles=2 on the three gap datasets from T-251.
#
# Usage:
#   Rscript dev/benchmarks/bench_drift_mpt.R

library(TreeSearch, lib.loc = ".agent-E")
library(TreeTools)
library(TreeDist)

DATASETS <- c("Wortley2006", "Zhu2013", "Geisler2001")
DRIFT_CONDITIONS <- c(0L, 2L)
SEEDS <- 1:3
BUDGETS <- c(30, 120)

# Use default preset parameters for everything except driftCycles.
# strategy = "none" bypasses auto-selection; explicit control overrides.
make_control <- function(drift_cycles) {
  SearchControl(
    tbrMaxHits = 1L,
    nniFirst = TRUE,
    sprFirst = FALSE,
    tabuSize = 100L,
    wagnerStarts = 3L,
    outerCycles = 1L,
    maxOuterResets = 2L,
    ratchetCycles = 12L,
    ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 0L,
    ratchetPerturbMaxMoves = 5L,
    adaptiveLevel = TRUE,
    driftCycles = drift_cycles,
    driftAfdLimit = 5L,
    driftRfdLimit = 0.15,
    xssRounds = 3L,
    xssPartitions = 4L,
    rssRounds = 1L,
    cssRounds = 0L,
    consensusStableReps = 3L,
    fuseInterval = 3L,
    fuseAcceptEqual = FALSE,
    poolMaxSize = 100L,
    enumTimeFraction = 0.1
  )
}

# Compute pairwise RF distances between trees, return summary stats
tree_diversity <- function(trees) {
  n <- length(trees)
  if (n < 2) return(list(mean_rf = NA, median_rf = NA, min_rf = NA, max_rf = NA))
  rf_mat <- as.matrix(RobinsonFoulds(trees))
  # Upper triangle only (exclude diagonal)
  rf_vals <- rf_mat[upper.tri(rf_mat)]
  list(
    mean_rf = mean(rf_vals),
    median_rf = median(rf_vals),
    min_rf = min(rf_vals),
    max_rf = max(rf_vals)
  )
}

results <- list()
row_i <- 0L

for (ds_name in DATASETS) {
  ds <- inapplicable.phyData[[ds_name]]
  n_tips <- length(ds)
  cat(sprintf("\n=== %s (%d tips) ===\n", ds_name, n_tips))

  for (budget in BUDGETS) {
    for (drift in DRIFT_CONDITIONS) {
      ctrl <- make_control(drift)
      for (seed in SEEDS) {
        row_i <- row_i + 1L
        cat(sprintf("  budget=%ds drift=%d seed=%d ... ", budget, drift, seed))
        t0 <- proc.time()

        res <- MaximizeParsimony(
          ds,
          maxSeconds = budget,
          strategy = "none",
          control = ctrl,
          verbosity = 0L,
          nThread = 1L
        )

        wall_s <- as.double((proc.time() - t0)[3])
        best_score <- attr(res, "score")
        n_trees <- length(res)
        n_topo <- attr(res, "n_topologies")
        n_reps <- attr(res, "replicates")
        timings <- attr(res, "timings")

        # Topological diversity (RF distances)
        div <- tree_diversity(res)

        cat(sprintf("score=%.0f trees=%d topo=%d reps=%d (%.1fs)\n",
                    best_score, n_trees, n_topo, n_reps, wall_s))

        results[[row_i]] <- data.frame(
          dataset = ds_name,
          n_tips = n_tips,
          budget_s = budget,
          drift_cycles = drift,
          seed = seed,
          best_score = best_score,
          n_trees = n_trees,
          n_topologies = n_topo,
          replicates = n_reps,
          wall_s = round(wall_s, 2),
          drift_ms = timings["drift_ms"],
          total_ms = sum(timings),
          drift_pct = round(100 * timings["drift_ms"] / sum(timings), 1),
          mean_rf = div$mean_rf,
          median_rf = div$median_rf,
          min_rf = div$min_rf,
          max_rf = div$max_rf,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

df <- do.call(rbind, results)
rownames(df) <- NULL

out_path <- "dev/benchmarks/results_drift_mpt.csv"
write.csv(df, out_path, row.names = FALSE)
cat(sprintf("\nResults written to %s\n", out_path))

# Quick summary table
cat("\n=== Summary by dataset × budget × drift ===\n")
agg <- aggregate(
  cbind(best_score, n_trees, n_topologies, replicates, mean_rf) ~ dataset + budget_s + drift_cycles,
  data = df, FUN = median
)
print(agg[order(agg$dataset, agg$budget_s, agg$drift_cycles), ])
