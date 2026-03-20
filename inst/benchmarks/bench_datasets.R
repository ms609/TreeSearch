# Benchmark dataset loading and scoring utilities
#
# Usage:
#   source("inst/benchmarks/bench_datasets.R")
#   datasets <- load_benchmark_datasets()
#   run_benchmark_suite(maxSeconds = 30, replicates = 5)

library(TreeSearch)
library(TreeTools)

# The 14 benchmark datasets, ordered by tip count
BENCHMARK_NAMES <- c(
 "Longrich2010",   # 20 tips, 3 states, 45% missing
  "Vinther2008",    # 23 tips, 4 states, 21% missing
  "Sansom2010",     # 23 tips, 4 states, 40% missing
  "DeAssis2011",    # 33 tips, 3 states, 21% inapp
  "Aria2015",       # 35 tips, 6 states, 13% missing
  "Wortley2006",    # 37 tips, 8 states, 31% missing
  "Griswold1999",   # 43 tips, 6 states, 6% missing
  "Schulze2007",    # 52 tips, 3 states, 17% inapp
  "Eklund2004",     # 54 tips, 6 states, 30% missing
  "Agnarsson2004",  # 62 tips, 7 states, 6% missing
  "Zanol2014",      # 74 tips, 9 states, 17% inapp
  "Zhu2013",        # 75 tips, 4 states, 43% missing
  "Giles2015",      # 78 tips, 4 states, 42% missing
  "Dikow2009"       # 88 tips, 9 states, 0.4% missing
)

#' Prepare raw data for C++ bridge from a phyDat object
#' @param dataset A phyDat object
#' @return List with contrast, tip_data, weight, levels
prepare_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    n_taxa = length(dataset)
  )
}

#' Load all 14 benchmark datasets
#' @return Named list of prepared datasets (ready for C++ bridge)
load_benchmark_datasets <- function() {
  datasets <- list()
  for (nm in BENCHMARK_NAMES) {
    ds <- TreeSearch::inapplicable.phyData[[nm]]
    if (is.null(ds)) {
      warning("Dataset ", nm, " not found in inapplicable.phyData")
      next
    }
    datasets[[nm]] <- prepare_ts_data(ds)
  }
  datasets
}

#' Characterize a benchmark dataset
#' @param ds Prepared dataset (from prepare_ts_data)
#' @return Data frame with one row of characteristics
characterize_dataset <- function(ds) {
  n_taxa <- ds$n_taxa
  n_patterns <- length(ds$weight)
  n_chars <- sum(ds$weight)
  lvls <- ds$levels
  contrast <- ds$contrast
  n_states <- ncol(contrast)
  inapp_idx <- which(lvls == "-")
  n_app_states <- n_states - length(inapp_idx)

  td <- ds$tip_data
  total_cells <- n_taxa * n_patterns

  n_inapp <- 0L
  n_missing <- 0L
  has_inapp <- length(inapp_idx) > 0
  for (i in seq_len(nrow(contrast))) {
    is_inapp <- has_inapp && contrast[i, inapp_idx] > 0.5
    cols_check <- setdiff(seq_len(n_states), inapp_idx)
    is_all <- all(contrast[i, cols_check] > 0.5)
    count <- sum(td == i)
    if (is_inapp && !is_all) n_inapp <- n_inapp + count
    if (is_all) n_missing <- n_missing + count
  }

  data.frame(
    n_taxa = n_taxa,
    n_chars = n_chars,
    n_patterns = n_patterns,
    pct_inapp = round(100 * n_inapp / total_cells, 1),
    n_app_states = n_app_states,
    pct_missing = round(100 * n_missing / total_cells, 1)
  )
}

#' Run a single benchmark: driven search on one dataset
#' @param name Dataset name (from BENCHMARK_NAMES)
#' @param maxSeconds Timeout in seconds
#' @param maxReplicates Maximum replicates
#' @param seed RNG seed
#' @param datasets Pre-loaded datasets (optional)
#' @return List with score, replicates, time, etc.
score_dataset <- function(name, maxSeconds = 10, maxReplicates = 20L,
                          seed = 42L, datasets = NULL) {
  if (is.null(datasets)) {
    ds <- prepare_ts_data(TreeSearch::inapplicable.phyData[[name]])
  } else {
    ds <- datasets[[name]]
  }
  if (is.null(ds)) stop("Dataset '", name, "' not found")

  set.seed(seed)
  t0 <- proc.time()
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = 5L,
    ratchetCycles = 5L,
    xssRounds = 1L,
    xssPartitions = 3L,
    fuseInterval = 5L,
    maxSeconds = maxSeconds,
    verbosity = 0L
  )
  elapsed <- (proc.time() - t0)[3]

  list(
    dataset = name,
    n_taxa = ds$n_taxa,
    best_score = result$best_score,
    replicates = result$replicates,
    pool_size = result$pool_size,
    hits_to_best = result$hits_to_best,
    timed_out = result$timed_out,
    elapsed = elapsed
  )
}

#' Run the full benchmark suite
#' @param maxSeconds Timeout per dataset
#' @param replicates Number of independent runs per dataset
#' @param seed Base seed (incremented per replicate)
#' @return Data frame with results
run_benchmark_suite <- function(maxSeconds = 30, replicates = 3L,
                                 seed = 42L) {
  datasets <- load_benchmark_datasets()
  results <- list()

  for (nm in names(datasets)) {
    for (rep in seq_len(replicates)) {
      cat(sprintf("[%s] rep %d/%d (timeout=%ds)...",
                  nm, rep, replicates, maxSeconds))
      res <- score_dataset(nm, maxSeconds = maxSeconds,
                            seed = seed + rep - 1L,
                            datasets = datasets)
      cat(sprintf(" score=%.0f reps=%d time=%.1fs\n",
                  res$best_score, res$replicates, res$elapsed))
      res$replicate <- rep
      results <- c(results, list(as.data.frame(res)))
    }
  }

  do.call(rbind, results)
}

#' Summarize benchmark results
#' @param results Data frame from run_benchmark_suite
#' @return Summary data frame (best score, median time, etc.)
summarize_benchmark <- function(results) {
  datasets <- unique(results$dataset)
  summaries <- list()

  for (nm in datasets) {
    sub <- results[results$dataset == nm, ]
    summaries <- c(summaries, list(data.frame(
      dataset = nm,
      n_taxa = sub$n_taxa[1],
      best_score = min(sub$best_score),
      median_score = median(sub$best_score),
      median_time = round(median(sub$elapsed), 2),
      median_reps = median(sub$replicates),
      stringsAsFactors = FALSE
    )))
  }

  do.call(rbind, summaries)
}
