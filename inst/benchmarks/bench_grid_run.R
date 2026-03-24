# Focused benchmark grid: no callback (workaround for segfault in progress_cb).
# Collects per-phase timings, wall-clock time, scores, convergence stats.

library(TreeSearch, lib.loc = if (dir.exists(".agent-a")) ".agent-a" else .libPaths())
library(TreeTools)

source("inst/benchmarks/bench_datasets.R")
source("inst/benchmarks/bench_framework.R")

# Simplified benchmark_run without callback
benchmark_run_nocb <- function(ds, strategy,
                                maxReplicates = 100L,
                                targetHits = NULL,
                                maxSeconds = 0,
                                seed = 42L) {
  if (is.null(targetHits)) {
    targetHits <- max(10L, ds$n_taxa %/% 5L)
  }

  args <- c(
    list(
      contrast = ds$contrast,
      tip_data = ds$tip_data,
      weight = ds$weight,
      levels = ds$levels,
      maxReplicates = as.integer(maxReplicates),
      targetHits = as.integer(targetHits),
      maxSeconds = as.double(maxSeconds),
      verbosity = 0L
    ),
    strategy
  )

  set.seed(seed)
  t0 <- proc.time()
  result <- do.call(TreeSearch:::ts_driven_search, args)
  wall_s <- as.double((proc.time() - t0)[3])

  list(
    best_score   = result$best_score,
    replicates   = result$replicates,
    hits_to_best = result$hits_to_best,
    pool_size    = result$pool_size,
    timed_out    = result$timed_out,
    wall_s       = wall_s,
    timings      = result$timings
  )
}

# Representative subset: small, medium, large datasets
GRID_DATASETS <- c(
  "Longrich2010",   # 20 tips
  "Vinther2008",    # 23 tips
  "Aria2015",       # 35 tips
  "Griswold1999",   # 43 tips
  "Agnarsson2004",  # 62 tips
  "Zhu2013",        # 75 tips
  "Giles2015",      # 78 tips
  "Dikow2009"       # 88 tips
)

run_grid <- function(dataset_names = GRID_DATASETS,
                     strategy_names = STRATEGY_NAMES,
                     replicates = 3L,
                     maxReplicates = 100L,
                     maxSeconds = 20,
                     base_seed = 7142L) {
  datasets <- load_benchmark_datasets()
  n_combos <- length(dataset_names) * length(strategy_names) * replicates
  cat(sprintf("Grid: %d datasets x %d strategies x %d reps = %d runs\n",
              length(dataset_names), length(strategy_names), replicates, n_combos))

  rows <- vector("list", n_combos)
  idx <- 0L

  for (ds_name in dataset_names) {
    ds <- datasets[[ds_name]]
    if (is.null(ds)) {
      warning("Skipping missing dataset: ", ds_name)
      next
    }
    for (strat_name in strategy_names) {
      strat <- get_strategy(strat_name)
      for (rep in seq_len(replicates)) {
        idx <- idx + 1L
        seed <- base_seed + (idx - 1L) * 7L

        cat(sprintf("[%3d/%d] %-15s x %-16s rep %d ...",
                    idx, n_combos, ds_name, strat_name, rep))

        res <- tryCatch(
          benchmark_run_nocb(ds, strat,
                             maxReplicates = maxReplicates,
                             targetHits = max(10L, ds$n_taxa %/% 5L),
                             maxSeconds = maxSeconds,
                             seed = seed),
          error = function(e) {
            cat(sprintf(" ERROR: %s\n", conditionMessage(e)))
            NULL
          }
        )

        if (is.null(res)) {
          rows[[idx]] <- data.frame(
            dataset = ds_name, strategy = strat_name, replicate = rep,
            seed = seed, n_taxa = ds$n_taxa,
            best_score = NA_real_, replicates = NA_integer_,
            hits_to_best = NA_integer_, pool_size = NA_integer_,
            timed_out = NA, wall_s = NA_real_,
            wagner_ms = NA_real_, tbr_ms = NA_real_,
            xss_ms = NA_real_, rss_ms = NA_real_, css_ms = NA_real_,
            ratchet_ms = NA_real_, drift_ms = NA_real_,
            final_tbr_ms = NA_real_, fuse_ms = NA_real_,
            stringsAsFactors = FALSE
          )
          next
        }

        cat(sprintf(" score=%.0f wall=%.1fs reps=%d %s\n",
                    res$best_score, res$wall_s, res$replicates,
                    if (res$timed_out) "[TIMEOUT]" else ""))

        rows[[idx]] <- data.frame(
          dataset = ds_name, strategy = strat_name, replicate = rep,
          seed = seed, n_taxa = ds$n_taxa,
          best_score = res$best_score, replicates = res$replicates,
          hits_to_best = res$hits_to_best, pool_size = res$pool_size,
          timed_out = res$timed_out, wall_s = res$wall_s,
          wagner_ms = res$timings[["wagner_ms"]],
          tbr_ms = res$timings[["tbr_ms"]],
          xss_ms = res$timings[["xss_ms"]],
          rss_ms = res$timings[["rss_ms"]],
          css_ms = res$timings[["css_ms"]],
          ratchet_ms = res$timings[["ratchet_ms"]],
          drift_ms = res$timings[["drift_ms"]],
          final_tbr_ms = res$timings[["final_tbr_ms"]],
          fuse_ms = res$timings[["fuse_ms"]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, rows[seq_len(idx)])
}

# Main
cat("Starting benchmark grid...\n\n")
results <- run_grid()
outfile <- "inst/benchmarks/results_grid.csv"
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nResults saved to %s (%d rows)\n", outfile, nrow(results)))
