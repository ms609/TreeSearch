# Brazeau-track baseline benchmark
#
# Benchmarks search strategies under Brazeau et al. (2019) NA-aware scoring
# on the fixed 20-matrix Brazeau sample. Compares strategy presets under
# both EW and IW(k=10), and includes a sample-size stability analysis.
#
# *** Run on Hamilton HPC, not locally. ***
#
# Usage:
#   source("dev/benchmarks/bench_brazeau_baseline.R")
#   results <- brazeau_benchmark_full()      # Hamilton only
#   results <- brazeau_benchmark_smoke()     # local quick test
#   stability <- sample_size_stability(results)
#
# Output: data frame with one row per dataset × weighting × strategy × seed.

library(TreeSearch)
library(TreeTools)

source("dev/benchmarks/bench_datasets.R")

# ---- Core benchmark runner ----

#' Run a single TreeSearch benchmark (Brazeau scoring)
#'
#' Uses MaximizeParsimony with default Brazeau scoring (no fitch_mode).
#' Returns score, timing, and replicate metadata.
#'
#' @param raw_dataset phyDat object (unmodified — Brazeau scoring applied)
#' @param timeout_s Timeout in seconds
#' @param seed RNG seed
#' @param concavity Concavity constant (Inf = EW, finite = IW)
#' @param strategy Strategy preset name
#' @param hits Target hits for convergence
#' @param reps Maximum replicates
#' @return Named list with score, timing, replicates, etc.
run_brazeau_search <- function(raw_dataset, timeout_s = 30, seed = 1,
                               concavity = Inf, strategy = "auto",
                               hits = 10L, reps = 50L) {
  set.seed(seed)
  t0 <- proc.time()
  result <- tryCatch(
    TreeSearch::MaximizeParsimony(
      raw_dataset,
      concavity = concavity,
      maxReplicates = as.integer(reps),
      targetHits = as.integer(hits),
      maxSeconds = as.double(timeout_s),
      strategy = strategy,
      verbosity = 0L,
      nThreads = 1L
    ),
    error = function(e) {
      warning("Search error: ", conditionMessage(e))
      structure(list(), class = "multiPhylo",
                score = NA_real_, replicates = NA_integer_,
                hits_to_best = NA_integer_, pool_size = NA_integer_)
    }
  )
  wall_s <- as.double((proc.time() - t0)[3])

  list(
    score = attr(result, "score"),
    n_trees = length(result),
    wall_s = wall_s,
    replicates = attr(result, "replicates"),
    hits = attr(result, "hits_to_best")
  )
}


# ---- Benchmark grid runner ----

#' Run the Brazeau benchmark grid
#'
#' @param dataset_keys Keys from MBANK_BRAZEAU_SAMPLE (or a subset)
#' @param weightings Named list: name -> concavity value
#' @param strategies Character vector of strategy presets
#' @param timeout_s Timeout per search
#' @param seeds Numeric vector of seeds
#' @param hits Target hits
#' @param reps Max replicates
#' @return data.frame with results
run_brazeau_grid <- function(
    dataset_keys = MBANK_BRAZEAU_SAMPLE,
    weightings = list(EW = Inf, IW10 = 10),
    strategies = c("default", "thorough"),
    timeout_s = 30,
    seeds = 1:5,
    hits = 10L,
    reps = 50L) {

  cat <- load_mbank_catalogue()
  datasets <- load_mbank_datasets(cat, dataset_keys, verbose = TRUE)

  # Load raw phyDat objects (need phyDat for MaximizeParsimony)
  raw_datasets <- list()
  for (k in names(datasets)) {
    row <- cat[cat$key == k, ]
    nex_path <- file.path(NEOTRANS_MATRICES_DIR, row$filename)
    raw_datasets[[k]] <- suppressWarnings(TreeTools::ReadAsPhyDat(nex_path))
  }

  n_combos <- length(dataset_keys) * length(weightings) * length(strategies) *
    length(seeds)
  cat(sprintf("\nBrazeau benchmark grid: %d combinations\n", n_combos))
  cat(sprintf("  Datasets: %d, Weightings: %s, Strategies: %s\n",
              length(dataset_keys),
              paste(names(weightings), collapse = "+"),
              paste(strategies, collapse = "+")))
  cat(sprintf("  Timeout: %ds, Seeds: %d\n\n", timeout_s, length(seeds)))

  rows <- list()
  idx <- 0L

  for (k in dataset_keys) {
    if (!k %in% names(raw_datasets)) next
    ds_info <- datasets[[k]]
    ds_raw <- raw_datasets[[k]]

    for (wt_name in names(weightings)) {
      conc <- weightings[[wt_name]]
      for (strat in strategies) {
        for (seed in seeds) {
          idx <- idx + 1L
          cat(sprintf("[%d/%d] %s %s %s seed=%d ... ",
                      idx, n_combos, k, wt_name, strat, seed))

          res <- run_brazeau_search(
            ds_raw, timeout_s = timeout_s, seed = seed,
            concavity = conc, strategy = strat,
            hits = hits, reps = reps
          )

          cat(sprintf("score=%.2f reps=%s time=%.1fs\n",
                      res$score,
                      ifelse(is.na(res$replicates), "NA",
                             as.character(res$replicates)),
                      res$wall_s))

          rows[[idx]] <- data.frame(
            dataset = k,
            n_taxa = ds_info$n_taxa,
            n_chars = sum(ds_info$weight),
            pct_inapp = cat[cat$key == k, "pct_inapp"],
            weighting = wt_name,
            concavity = conc,
            strategy = strat,
            timeout_s = timeout_s,
            seed = seed,
            score = res$score,
            n_trees = res$n_trees,
            wall_s = round(res$wall_s, 3),
            replicates = res$replicates,
            hits = res$hits,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  do.call(rbind, rows)
}


# ---- Convenience wrappers ----

#' Smoke test: 3 datasets, EW only, 5s
brazeau_benchmark_smoke <- function() {
  run_brazeau_grid(
    dataset_keys = MBANK_BRAZEAU_SAMPLE[c(1, 7, 14)],
    weightings = list(EW = Inf),
    strategies = c("default"),
    timeout_s = 5,
    seeds = 1:2,
    hits = 3L,
    reps = 10L
  )
}

#' Full benchmark: all 20 datasets, EW + IW(k=10), default + thorough
brazeau_benchmark_full <- function() {
  run_brazeau_grid(
    dataset_keys = MBANK_BRAZEAU_SAMPLE,
    weightings = list(EW = Inf, IW10 = 10),
    strategies = c("default", "thorough"),
    timeout_s = 30,
    seeds = 1:5,
    hits = 10L,
    reps = 50L
  )
}


# ---- Sample-size stability analysis ----

#' Assess whether 20 matrices are sufficient for stable strategy ranking
#'
#' Bootstraps subsamples of increasing size from benchmark results,
#' ranks strategies within each subsample, and computes rank concordance
#' with the full-sample ranking.
#'
#' @param results data.frame from run_brazeau_grid (single weighting)
#' @param n_boot Number of bootstrap resamples per k
#' @param ks Subsample sizes to test
#' @return data.frame with columns: k, median_W, q25_W, q75_W
sample_size_stability <- function(results,
                                  n_boot = 500L,
                                  ks = c(5, 8, 10, 12, 14, 16, 18, 20)) {

  # Compute per-dataset × strategy mean score
  agg <- aggregate(score ~ dataset + strategy, data = results, FUN = mean)
  datasets <- unique(agg$dataset)
  strategies <- unique(agg$strategy)

  if (length(strategies) < 2) {
    stop("Need at least 2 strategies for ranking stability analysis")
  }

  # Compute score gap: strategy_score - best_strategy_score, per dataset
  best_per_ds <- aggregate(score ~ dataset, data = agg, FUN = min)
  names(best_per_ds)[2] <- "best_score"
  agg <- merge(agg, best_per_ds, by = "dataset")
  agg$gap <- agg$score - agg$best_score

  # Full-sample ranking (lower gap = better)
  full_rank <- aggregate(gap ~ strategy, data = agg, FUN = mean)
  full_rank <- full_rank[order(full_rank$gap), ]
  full_ranking <- seq_len(nrow(full_rank))
  names(full_ranking) <- full_rank$strategy

  # Bootstrap at each k
  stability <- data.frame(k = integer(), median_W = numeric(),
                          q25_W = numeric(), q75_W = numeric())

  for (k in ks) {
    if (k > length(datasets)) next

    W_values <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      sampled_ds <- sample(datasets, k, replace = FALSE)
      sub <- agg[agg$dataset %in% sampled_ds, ]
      sub_rank <- aggregate(gap ~ strategy, data = sub, FUN = mean)
      sub_rank <- sub_rank[order(sub_rank$gap), ]

      # Kendall's tau-b between full and subsample rankings
      sub_ranking <- match(names(full_ranking), sub_rank$strategy)
      tau <- cor(full_ranking, sub_ranking, method = "kendall",
                 use = "complete.obs")
      W_values[b] <- (tau + 1) / 2  # normalize to [0, 1]
    }

    stability <- rbind(stability, data.frame(
      k = k,
      median_W = median(W_values, na.rm = TRUE),
      q25_W = quantile(W_values, 0.25, na.rm = TRUE),
      q75_W = quantile(W_values, 0.75, na.rm = TRUE)
    ))
  }

  stability
}


# ---- Results summary ----

#' Summarize Brazeau benchmark results
#'
#' @param results data.frame from run_brazeau_grid
#' @return Summary data.frame with per-dataset × weighting × strategy stats
summarize_brazeau <- function(results) {
  results |>
    aggregate(cbind(score, wall_s, replicates) ~
                dataset + n_taxa + weighting + strategy,
              data = _, FUN = function(x) {
                c(best = min(x, na.rm = TRUE),
                  median = median(x, na.rm = TRUE))
              }) |>
    do.call(data.frame, args = _)
}

#' Save results to CSV
save_brazeau_results <- function(results,
                                 file = sprintf(
                                   "dev/benchmarks/brazeau_baseline_%s.csv",
                                   format(Sys.time(), "%Y%m%d_%H%M"))) {
  write.csv(results, file, row.names = FALSE)
  cat("Results saved to", file, "\n")
  invisible(file)
}
