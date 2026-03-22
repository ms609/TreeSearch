# Phase 6D: Benchmarking framework
#
# Runs dataset x strategy x N replicates and records:
#   - Best score found
#   - Total wall-clock time
#   - Time to best score (via progress callback)
#   - Number of replicates to convergence
#   - Per-phase timing breakdown
#
# Usage:
#   source("inst/benchmarks/bench_framework.R")
#   results <- run_benchmark_grid()
#   summary <- summarize_grid(results)

library(TreeSearch)
library(TreeTools)

source("inst/benchmarks/bench_datasets.R")

# ---- Strategy presets (formalized from strategies.md, T-003) ----

STRATEGY_NAMES <- c("sprint", "default", "thorough",
                    "ratchet_heavy", "sectorial_heavy", "drift_heavy")

get_strategy <- function(name = STRATEGY_NAMES) {
  name <- match.arg(name)
  strategies <- list(
    sprint = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 0L,
      ratchetCycles = 3L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 0L, driftAfdLimit = 3L, driftRfdLimit = 0.1,
      xssRounds = 1L, xssPartitions = 4L, rssRounds = 0L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 5L, fuseAcceptEqual = FALSE
    ),
    default = list(
      wagnerStarts = 3L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 10L, ratchetPerturbProb = 0.25,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 5L,
      ratchetAdaptive = FALSE,
      driftCycles = 4L, driftAfdLimit = 5L, driftRfdLimit = 0.15,
      xssRounds = 3L, xssPartitions = 4L, rssRounds = 1L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 3L, fuseAcceptEqual = FALSE,
      sprFirst = TRUE, adaptiveLevel = TRUE, consensusStableReps = 3L
    ),
    thorough = list(
      wagnerStarts = 3L, tbrMaxHits = 3L, tabuSize = 200L,
      ratchetCycles = 20L, ratchetPerturbProb = 0.25,
      ratchetPerturbMode = 2L, ratchetPerturbMaxMoves = 5L,
      ratchetAdaptive = TRUE,
      driftCycles = 12L, driftAfdLimit = 5L, driftRfdLimit = 0.15,
      xssRounds = 5L, xssPartitions = 6L, rssRounds = 3L,
      cssRounds = 2L, cssPartitions = 6L,
      sectorMinSize = 6L, sectorMaxSize = 80L,
      fuseInterval = 2L, fuseAcceptEqual = TRUE
    ),
    ratchet_heavy = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 30L, ratchetPerturbProb = 0.30,
      ratchetPerturbMode = 2L, ratchetPerturbMaxMoves = 5L,
      ratchetAdaptive = TRUE,
      driftCycles = 2L, driftAfdLimit = 3L, driftRfdLimit = 0.1,
      xssRounds = 1L, xssPartitions = 4L, rssRounds = 0L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 3L, fuseAcceptEqual = FALSE
    ),
    sectorial_heavy = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 5L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 3L, driftAfdLimit = 3L, driftRfdLimit = 0.1,
      xssRounds = 8L, xssPartitions = 6L, rssRounds = 4L,
      cssRounds = 3L, cssPartitions = 6L,
      sectorMinSize = 6L, sectorMaxSize = 80L,
      fuseInterval = 2L, fuseAcceptEqual = TRUE
    ),
    drift_heavy = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 5L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 20L, driftAfdLimit = 5L, driftRfdLimit = 0.2,
      xssRounds = 2L, xssPartitions = 4L, rssRounds = 1L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 3L, fuseAcceptEqual = TRUE
    )
  )
  strategies[[name]]
}

# ---- Best-known EW scores (from datasets.md, T-002) ----

BEST_KNOWN_EW <- c(
  Longrich2010 = 131, Vinther2008 = 79, Sansom2010 = 189,
  DeAssis2011 = 64, Aria2015 = 145, Wortley2006 = 496,
  Griswold1999 = 409, Schulze2007 = 167, Eklund2004 = 445,
  Agnarsson2004 = 778, Zanol2014 = 1338, Zhu2013 = 649,
  Giles2015 = 720, Dikow2009 = 1614
)

# ---- Core benchmark function ----

#' Run one driven search and record performance metrics.
#'
#' Calls ts_driven_search directly with the given strategy parameters.
#' Uses a progress callback to record the wall-clock time at which the
#' best score was first found ("time to best").
#'
#' @param ds Prepared dataset (from prepare_ts_data).
#' @param strategy Named list of strategy parameters (from get_strategy).
#' @param maxReplicates Hard replicate cap.
#' @param targetHits Convergence criterion (hits to best score).
#' @param maxSeconds Wall-clock timeout (0 = no timeout).
#' @param seed RNG seed.
#' @return Named list with score, timing, and convergence metrics.
benchmark_run <- function(ds, strategy,
                          maxReplicates = 100L,
                          targetHits = NULL,
                          maxSeconds = 0,
                          seed = 42L) {
  if (is.null(targetHits)) {
    targetHits <- max(10L, ds$n_taxa %/% 5L)
  }

  # Progress-callback state: track when best score first appeared
  cb_env <- new.env(parent = emptyenv())
  cb_env$best <- Inf
  cb_env$time_to_best <- NA_real_
  cb_env$trace <- list()

  progress_cb <- function(info) {
    if (is.finite(info$best_score) && info$best_score < cb_env$best) {
      cb_env$best <- info$best_score
      cb_env$time_to_best <- info$elapsed
    }
    cb_env$trace[[length(cb_env$trace) + 1L]] <- list(
      replicate = info$replicate,
      elapsed = info$elapsed,
      best_score = info$best_score,
      hits = info$hits_to_best,
      phase = info$phase
    )
  }

  # Build the full argument list for ts_driven_search.
  # verbosity >= 1 required for the C++ engine to invoke the callback.
  args <- c(
    list(
      contrast = ds$contrast,
      tip_data = ds$tip_data,
      weight = ds$weight,
      levels = ds$levels,
      maxReplicates = as.integer(maxReplicates),
      targetHits = as.integer(targetHits),
      maxSeconds = as.double(maxSeconds),
      verbosity = 1L,
      progressCallback = progress_cb
    ),
    strategy
  )

  set.seed(seed)
  t0 <- proc.time()
  result <- do.call(TreeSearch:::ts_driven_search, args)
  wall_s <- as.double((proc.time() - t0)[3])

  list(
    best_score     = result$best_score,
    replicates     = result$replicates,
    hits_to_best   = result$hits_to_best,
    pool_size      = result$pool_size,
    timed_out      = result$timed_out,
    wall_s         = wall_s,
    time_to_best_s = cb_env$time_to_best,
    timings        = result$timings,
    trace          = cb_env$trace
  )
}

# ---- Grid runner ----

#' Run the full dataset x strategy x replicate benchmark grid.
#'
#' @param dataset_names Character vector of dataset names.
#' @param strategy_names Character vector of strategy preset names.
#' @param replicates Number of independent runs per combination.
#' @param maxReplicates Replicate cap per run.
#' @param targetHits Convergence hits (NULL = auto).
#' @param maxSeconds Timeout per run (0 = no timeout).
#' @param base_seed Seed for first replicate; incremented per replicate.
#' @return A data.frame with one row per dataset x strategy x replicate.
run_benchmark_grid <- function(
    dataset_names = BENCHMARK_NAMES,
    strategy_names = STRATEGY_NAMES,
    replicates = 5L,
    maxReplicates = 100L,
    targetHits = NULL,
    maxSeconds = 30,
    base_seed = 42L
) {
  datasets <- load_benchmark_datasets()
  n_combos <- length(dataset_names) * length(strategy_names) * replicates
  cat(sprintf("Benchmark grid: %d datasets x %d strategies x %d reps = %d runs\n",
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
        seed <- base_seed + rep - 1L

        cat(sprintf("[%3d/%d] %s x %s rep %d ...",
                    idx, n_combos, ds_name, strat_name, rep))

        res <- tryCatch(
          benchmark_run(ds, strat,
                        maxReplicates = maxReplicates,
                        targetHits = targetHits,
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
            time_to_best_s = NA_real_,
            wagner_ms = NA_real_, tbr_ms = NA_real_,
            xss_ms = NA_real_, rss_ms = NA_real_, css_ms = NA_real_,
            ratchet_ms = NA_real_, drift_ms = NA_real_,
            final_tbr_ms = NA_real_, fuse_ms = NA_real_,
            stringsAsFactors = FALSE
          )
          next
        }

        cat(sprintf(" score=%.0f wall=%.1fs ttb=%.1fs reps=%d\n",
                    res$best_score, res$wall_s,
                    if (is.na(res$time_to_best_s)) -1 else res$time_to_best_s,
                    res$replicates))

        rows[[idx]] <- data.frame(
          dataset = ds_name,
          strategy = strat_name,
          replicate = rep,
          seed = seed,
          n_taxa = ds$n_taxa,
          best_score = res$best_score,
          replicates = res$replicates,
          hits_to_best = res$hits_to_best,
          pool_size = res$pool_size,
          timed_out = res$timed_out,
          wall_s = res$wall_s,
          time_to_best_s = res$time_to_best_s,
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

# ---- Summarization ----

#' Summarize benchmark grid results per dataset x strategy.
#'
#' Computes: best score, median score, convergence rate (fraction that
#' hit targetHits before timeout), median wall time, median time-to-best,
#' and per-phase time medians.
#'
#' @param results Data frame from run_benchmark_grid.
#' @param best_known Named numeric vector of best-known EW scores.
#' @return Data frame with one row per dataset x strategy.
summarize_grid <- function(results, best_known = BEST_KNOWN_EW) {
  combos <- unique(results[, c("dataset", "strategy")])
  out <- vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    ds_name <- combos$dataset[i]
    st_name <- combos$strategy[i]
    sub <- results[results$dataset == ds_name & results$strategy == st_name, ]
    sub <- sub[!is.na(sub$best_score), , drop = FALSE]

    if (nrow(sub) == 0) next

    bk <- if (ds_name %in% names(best_known)) best_known[[ds_name]] else NA_real_

    # How many runs found the best-known score?
    found_optimal <- if (is.na(bk)) NA_real_ else mean(sub$best_score <= bk)

    total_phase_ms <- sub$wagner_ms + sub$tbr_ms + sub$xss_ms + sub$rss_ms +
      sub$css_ms + sub$ratchet_ms + sub$drift_ms + sub$final_tbr_ms +
      sub$fuse_ms

    out[[i]] <- data.frame(
      dataset = ds_name,
      strategy = st_name,
      n_taxa = sub$n_taxa[1],
      n_runs = nrow(sub),
      best_score = min(sub$best_score),
      median_score = median(sub$best_score),
      best_known = if (is.na(bk)) NA_real_ else bk,
      pct_found_optimal = round(100 * found_optimal, 1),
      converge_rate = round(100 * mean(!sub$timed_out), 1),
      median_wall_s = round(median(sub$wall_s), 3),
      median_ttb_s = round(median(sub$time_to_best_s, na.rm = TRUE), 3),
      median_reps = median(sub$replicates),
      median_hits = median(sub$hits_to_best),
      # Phase fraction (median % of total C++ time)
      pct_wagner = round(100 * median(sub$wagner_ms / total_phase_ms,
                                       na.rm = TRUE), 1),
      pct_tbr = round(100 * median(sub$tbr_ms / total_phase_ms,
                                    na.rm = TRUE), 1),
      pct_xss = round(100 * median(sub$xss_ms / total_phase_ms,
                                    na.rm = TRUE), 1),
      pct_rss = round(100 * median(sub$rss_ms / total_phase_ms,
                                    na.rm = TRUE), 1),
      pct_css = round(100 * median(sub$css_ms / total_phase_ms,
                                    na.rm = TRUE), 1),
      pct_ratchet = round(100 * median(sub$ratchet_ms / total_phase_ms,
                                        na.rm = TRUE), 1),
      pct_drift = round(100 * median(sub$drift_ms / total_phase_ms,
                                      na.rm = TRUE), 1),
      pct_fuse = round(100 * median(sub$fuse_ms / total_phase_ms,
                                     na.rm = TRUE), 1),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

# ---- Persistence helpers ----

#' Save benchmark results to CSV.
save_results <- function(results,
                         file = sprintf("inst/benchmarks/results_%s.csv",
                                        format(Sys.time(), "%Y%m%d_%H%M"))) {
  write.csv(results, file, row.names = FALSE)
  cat("Results saved to", file, "\n")
  invisible(file)
}

#' Load benchmark results from CSV.
load_results <- function(file) {
  read.csv(file, stringsAsFactors = FALSE)
}

# ---- Quick-start convenience wrappers ----

#' Run a small smoke test: 2 datasets x 2 strategies x 2 reps, 5s timeout.
benchmark_smoke <- function() {
  run_benchmark_grid(
    dataset_names = c("Vinther2008", "Agnarsson2004"),
    strategy_names = c("sprint", "default"),
    replicates = 2L,
    maxReplicates = 20L,
    maxSeconds = 5,
    base_seed = 42L
  )
}

#' Run the full production benchmark (all 14 datasets x 6 strategies).
#'
#' Warning: this takes a long time. At 30s timeout per run with 5 reps:
#' 14 x 6 x 5 = 420 runs x 30s = ~3.5 hours worst case.
benchmark_full <- function(maxSeconds = 30, replicates = 5L) {
  run_benchmark_grid(
    maxReplicates = 100L,
    maxSeconds = maxSeconds,
    replicates = replicates,
    base_seed = 42L
  )
}
