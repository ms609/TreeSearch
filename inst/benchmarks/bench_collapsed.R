#!/usr/bin/env Rscript
# Benchmark collapsed-tree optimization: skip counts, wall time, score equivalence
#
# Usage: Rscript inst/benchmarks/bench_collapsed.R <lib_path>
#
# Runs each dataset 3 times with fixed seeds and reports:
#   - Skip counts (via ts_tbr_search on near-optimal tree)
#   - Driven search wall time and score
#   - Per-phase timing breakdown

args <- commandArgs(trailingOnly = TRUE)
lib_path <- if (length(args) >= 1) args[1] else ".agent-a"

library(TreeSearch, lib.loc = lib_path)
library(TreeTools)

# --- Datasets ---
datasets <- c("Vinther2008", "Agnarsson2004", "Zhu2013", "Dikow2009")

prepare <- function(name) {
  ds <- TreeSearch::inapplicable.phyData[[name]]
  at <- attributes(ds)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(ds, use.names = FALSE),
                      nrow = length(ds), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    n_taxa = length(ds)
  )
}

# --- Part 1: Skip count measurement via TBR on near-optimal trees ---
cat("=== Part 1: Collapsed-flag skip counts (TBR) ===\n\n")

for (nm in datasets) {
  d <- prepare(nm)
  n_tip <- d$n_taxa
  n_internal <- n_tip - 1L
  total_clips <- n_tip + n_internal - 1L  # all nodes except root

  # Build a near-optimal tree via short driven search
  set.seed(7391)
  quick <- TreeSearch:::ts_driven_search(
    d$contrast, d$tip_data, d$weight, d$levels,
    maxReplicates = 3L, targetHits = 2L,
    ratchetCycles = 3L, driftCycles = 1L,
    xssRounds = 1L, xssPartitions = 3L,
    rssRounds = 0L, cssRounds = 0L,
    fuseInterval = 3L, maxSeconds = 30,
    verbosity = 0L, nThreads = 1L
  )

  # Run TBR from that tree (converged = already at local optimum)
  edge <- quick$trees[[1]]
  tbr_res <- TreeSearch:::ts_tbr_search(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = 10L, acceptEqual = FALSE
  )

  pct_skip <- round(100 * tbr_res$n_zero_skipped /
    (tbr_res$n_zero_skipped + tbr_res$n_evaluated), 1)

  cat(sprintf("%-15s  tips=%d  score=%.0f  evaluated=%d  skipped=%d  skip%%=%.1f%%\n",
              nm, n_tip, tbr_res$score,
              tbr_res$n_evaluated, tbr_res$n_zero_skipped, pct_skip))
}

# --- Part 2: Driven search wall time & scores ---
cat("\n=== Part 2: Driven search (3 seeds × 4 datasets) ===\n\n")

seeds <- c(2847L, 5192L, 8634L)
results <- list()

for (nm in datasets) {
  d <- prepare(nm)

  for (s in seeds) {
    set.seed(s)
    t0 <- proc.time()
    res <- TreeSearch:::ts_driven_search(
      d$contrast, d$tip_data, d$weight, d$levels,
      maxReplicates = 10L, targetHits = 5L,
      ratchetCycles = 5L, driftCycles = 2L,
      xssRounds = 3L, xssPartitions = 4L,
      rssRounds = 1L, cssRounds = 0L,
      fuseInterval = 3L, maxSeconds = 60,
      verbosity = 0L, nThreads = 1L
    )
    elapsed <- (proc.time() - t0)[3]

    tim <- res$timings
    row <- data.frame(
      dataset = nm,
      seed = s,
      score = res$best_score,
      reps = res$replicates,
      pool = res$pool_size,
      wall_s = round(elapsed, 2),
      tbr_ms = round(tim[["tbr_ms"]], 0),
      ratchet_ms = round(tim[["ratchet_ms"]], 0),
      drift_ms = round(tim[["drift_ms"]], 0),
      xss_ms = round(tim[["xss_ms"]], 0),
      rss_ms = round(tim[["rss_ms"]], 0),
      fuse_ms = round(tim[["fuse_ms"]], 0),
      final_tbr_ms = round(tim[["final_tbr_ms"]], 0),
      stringsAsFactors = FALSE
    )
    results <- c(results, list(row))
    cat(sprintf("  %-15s seed=%d  score=%.0f  reps=%d  wall=%.2fs\n",
                nm, s, res$best_score, res$replicates, elapsed))
  }
}

results_df <- do.call(rbind, results)

cat("\n=== Summary by dataset ===\n\n")
for (nm in datasets) {
  sub <- results_df[results_df$dataset == nm, ]
  cat(sprintf("%-15s  best=%.0f  median_wall=%.2fs  median_tbr_ms=%.0f  median_ratchet_ms=%.0f  median_drift_ms=%.0f\n",
              nm,
              min(sub$score),
              median(sub$wall_s),
              median(sub$tbr_ms),
              median(sub$ratchet_ms),
              median(sub$drift_ms)))
}

cat("\nDone.\n")
