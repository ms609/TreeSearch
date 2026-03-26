# Warm-start benchmark: measure ratchet/drift escape effectiveness
#
# Seeds search with a pre-computed TBR-optimal tree to isolate
# perturbation quality from initial descent quality.
#
# Usage:
#   source("dev/benchmarks/bench_framework.R")
#   source("dev/benchmarks/bench_warmstart.R")
#   ws <- warmstart_benchmark("Agnarsson2004", replicates = 20)
#   warmstart_summary(ws)

library(TreeSearch)
library(TreeTools)

source("dev/benchmarks/bench_datasets.R")

#' Compute a TBR-optimal tree via a short sprint search.
#'
#' Runs a fast search (sprint strategy, 1 replicate) to produce a local
#' optimum. This tree serves as the warm-start seed for escape benchmarks.
#'
#' @param ds Prepared dataset (from prepare_ts_data).
#' @param seed RNG seed for the sprint search.
#' @return Named list with `edge` (edge matrix), `score` (optimum score),
#'   and `tree` (phylo object) for inspection.
compute_warmstart_tree <- function(ds, seed = 7381L) {
  set.seed(seed)
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 1L,
    targetHits = 1L,
    ratchetCycles = 0L,
    driftCycles = 0L,
    xssRounds = 0L,
    rssRounds = 0L,
    cssRounds = 0L,
    nniPerturbCycles = 0L,
    maxSeconds = 0,
    verbosity = 0L
  )
  if (result$pool_size == 0) stop("Sprint search produced no trees")

  edge_mat <- result$trees[[1]]
  list(
    edge = edge_mat,
    score = result$best_score
  )
}

#' Run a single warm-started replicate.
#'
#' Passes the pre-computed tree via `startEdge`, runs 1 replicate with
#' the given strategy. Since the starting tree is already TBR-optimal,
#' the initial TBR phase converges immediately; only ratchet/drift/XSS
#' perturbations can improve the score.
#'
#' @param ds Prepared dataset.
#' @param start_edge Edge matrix from compute_warmstart_tree().
#' @param strategy Named list of strategy params (from get_strategy).
#' @param seed RNG seed for this replicate.
#' @param maxSeconds Timeout.
#' @return Named list with metrics.
warmstart_run <- function(ds, start_edge, strategy,
                          seed = 42L, maxSeconds = 30) {

  # Track when the score improves from the warm-start baseline
  cb_env <- new.env(parent = emptyenv())
  cb_env$best <- Inf
  cb_env$time_to_improvement <- NA_real_
  cb_env$trace <- list()

  progress_cb <- function(info) {
    if (is.finite(info$best_score) && info$best_score < cb_env$best) {
      cb_env$best <- info$best_score
      cb_env$time_to_improvement <- info$elapsed
    }
    cb_env$trace[[length(cb_env$trace) + 1L]] <- list(
      replicate = info$replicate,
      elapsed = info$elapsed,
      best_score = info$best_score,
      phase = info$phase
    )
  }

  args <- c(
    list(
      contrast = ds$contrast,
      tip_data = ds$tip_data,
      weight = ds$weight,
      levels = ds$levels,
      maxReplicates = 1L,
      targetHits = 1L,
      maxSeconds = as.double(maxSeconds),
      verbosity = 1L,
      startEdge = start_edge,
      progressCallback = progress_cb
    ),
    strategy
  )

  set.seed(seed)
  t0 <- proc.time()
  result <- do.call(TreeSearch:::ts_driven_search, args)
  wall_s <- as.double((proc.time() - t0)[3])

  list(
    best_score = result$best_score,
    wall_s = wall_s,
    time_to_improvement_s = cb_env$time_to_improvement,
    timed_out = result$timed_out,
    timings = result$timings,
    trace = cb_env$trace
  )
}

#' Run warm-start escape benchmark for one dataset.
#'
#' First computes a TBR-local-optimum via sprint, then runs multiple
#' warm-started replicates with varying seeds and strategies.
#'
#' @param ds_name Dataset name (from BENCHMARK_NAMES or LARGE_BENCHMARK_NAMES).
#' @param strategy_names Strategies to test.
#' @param replicates Independent warm-started runs per strategy.
#' @param maxSeconds Timeout per run.
#' @param warmstart_seed Seed for the initial sprint search.
#' @param base_seed Base seed for warm-started replicates.
#' @return Data frame with one row per strategy x replicate.
warmstart_benchmark <- function(
    ds_name,
    strategy_names = c("default", "thorough"),
    replicates = 10L,
    maxSeconds = 30,
    warmstart_seed = 7381L,
    base_seed = 42L
) {
  all_ds <- load_all_benchmark_datasets()
  ds <- all_ds[[ds_name]]
  if (is.null(ds)) stop("Dataset '", ds_name, "' not found")

  cat(sprintf("Computing warm-start tree for %s (%d tips)...\n",
              ds_name, ds$n_taxa))
  ws <- compute_warmstart_tree(ds, seed = warmstart_seed)
  cat(sprintf("Warm-start score: %.5g\n\n", ws$score))

  rows <- list()
  for (strat_name in strategy_names) {
    strat <- get_strategy(strat_name)
    for (rep in seq_len(replicates)) {
      seed <- base_seed + rep - 1L
      cat(sprintf("[%s rep %d/%d] ...", strat_name, rep, replicates))

      res <- tryCatch(
        warmstart_run(ds, ws$edge, strat, seed = seed,
                      maxSeconds = maxSeconds),
        error = function(e) {
          cat(sprintf(" ERROR: %s\n", conditionMessage(e)))
          NULL
        }
      )

      if (is.null(res)) {
        rows <- c(rows, list(data.frame(
          dataset = ds_name, n_taxa = ds$n_taxa,
          strategy = strat_name, replicate = rep, seed = seed,
          warmstart_score = ws$score,
          best_score = NA_real_, improvement = NA_real_,
          wall_s = NA_real_, time_to_improvement_s = NA_real_,
          timed_out = NA,
          stringsAsFactors = FALSE
        )))
        next
      }

      improvement <- ws$score - res$best_score
      cat(sprintf(" score=%.5g improvement=%.5g time=%.1fs\n",
                  res$best_score, improvement, res$wall_s))

      rows <- c(rows, list(data.frame(
        dataset = ds_name, n_taxa = ds$n_taxa,
        strategy = strat_name, replicate = rep, seed = seed,
        warmstart_score = ws$score,
        best_score = res$best_score,
        improvement = improvement,
        wall_s = res$wall_s,
        time_to_improvement_s = res$time_to_improvement_s,
        timed_out = res$timed_out,
        stringsAsFactors = FALSE
      )))
    }
  }

  do.call(rbind, rows)
}

#' Summarize warm-start benchmark results.
#'
#' @param results Data frame from warmstart_benchmark.
#' @return Summary per strategy: median improvement, escape rate, timing.
warmstart_summary <- function(results) {
  strats <- unique(results$strategy)
  summaries <- list()
  for (st in strats) {
    sub <- results[results$strategy == st & !is.na(results$best_score), ]
    if (nrow(sub) == 0) next
    escaped <- sub$improvement > 0
    summaries <- c(summaries, list(data.frame(
      strategy = st,
      n_runs = nrow(sub),
      warmstart_score = sub$warmstart_score[1],
      best_found = min(sub$best_score),
      median_score = median(sub$best_score),
      median_improvement = median(sub$improvement),
      escape_rate = round(100 * mean(escaped), 1),
      median_wall_s = round(median(sub$wall_s), 2),
      median_tti_s = round(median(sub$time_to_improvement_s, na.rm = TRUE), 2),
      stringsAsFactors = FALSE
    )))
  }
  do.call(rbind, summaries)
}
