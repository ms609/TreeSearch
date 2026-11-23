#!/usr/bin/env Rscript
# Performance benchmark for C++ driven search engine.
# Generates baseline timings for Phase 3 optimization work.
#
# Usage:
#   Rscript tests/benchmark/bench-driven-search.R [output_file]
#
# Output: CSV file with per-dataset timing breakdown.

.libPaths(c(".agent-e", .libPaths()))
library(TreeSearch)

# Representative datasets: small (20-35 tips), medium (50-68 tips), large (74-88 tips)
BENCH_DATASETS <- list(
  # Small
  small_na   = list(name = "Vinther2008",  tips = 23, chars = 50),
  small_na2  = list(name = "Longrich2010", tips = 20, chars = 80),
  small_na3  = list(name = "Sano2011",     tips = 36, chars = 52),
  # Medium
  med_na     = list(name = "Eklund2004",   tips = 54, chars = 131),
  med_na2    = list(name = "Wilson2003",    tips = 61, chars = 161),
  med_na3    = list(name = "Conrad2008",    tips = 64, chars = 360),
  # Large
  large_na   = list(name = "Zanol2014",    tips = 74, chars = 210),
  large_na2  = list(name = "Zhu2013",      tips = 75, chars = 253),
  large_na3  = list(name = "Dikow2009",    tips = 88, chars = 204)
)

# Benchmark parameters
N_REPS <- 3  # replicates per configuration (for timing stability)
MAX_SECONDS <- 30  # timeout per run
TARGET_HITS <- 3

run_benchmark <- function(dataset_name, mode = "EW", reps = N_REPS,
                          max_seconds = MAX_SECONDS) {
  dat <- inapplicable.phyData[[dataset_name]]
  if (is.null(dat)) stop("Dataset not found: ", dataset_name)

  concavity_val <- switch(mode,
    EW = Inf,
    IW3 = 3,
    IW10 = 10,
    stop("Unknown mode: ", mode)
  )

  times <- numeric(reps)
  scores <- numeric(reps)
  n_replicates <- integer(reps)
  pool_sizes <- integer(reps)
  timed_out <- logical(reps)

  for (i in seq_len(reps)) {
    set.seed(7291 + i)
    t0 <- proc.time()["elapsed"]
    result <- MaximizeParsimony(
      dat,
      concavity = concavity_val,
      maxReplicates = 50L,
      targetHits = TARGET_HITS,
      verbosity = 0L
    )
    t1 <- proc.time()["elapsed"]

    times[i] <- t1 - t0
    scores[i] <- attr(result, "score")
    n_replicates[i] <- attr(result, "replicates")
    pool_sizes[i] <- length(result)
    timed_out[i] <- isTRUE(attr(result, "timed_out"))
  }

  data.frame(
    dataset = dataset_name,
    n_tip = length(dat),
    n_char = attr(dat, "nr"),
    mode = mode,
    median_time = median(times),
    min_time = min(times),
    max_time = max(times),
    best_score = min(scores),
    median_replicates = median(n_replicates),
    median_pool = median(pool_sizes),
    any_timeout = any(timed_out),
    stringsAsFactors = FALSE
  )
}

# Run benchmarks
cat("TreeSearch driven search benchmark\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("TreeSearch version: %s\n", packageVersion("TreeSearch")))
cat("---\n\n")

results <- list()
idx <- 0

for (ds_info in BENCH_DATASETS) {
  for (mode in c("EW", "IW10")) {
    idx <- idx + 1
    cat(sprintf("[%2d/%2d] %-15s %-4s ... ",
                idx, length(BENCH_DATASETS) * 2, ds_info$name, mode))

    res <- tryCatch(
      run_benchmark(ds_info$name, mode),
      error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        NULL
      }
    )

    if (!is.null(res)) {
      cat(sprintf("%.2fs (score=%.1f, reps=%d)\n",
                  res$median_time, res$best_score, res$median_replicates))
      results[[idx]] <- res
    }
  }
}

bench_df <- do.call(rbind, results)

# Print summary table
cat("\n=== BENCHMARK RESULTS ===\n\n")
print(bench_df[, c("dataset", "n_tip", "n_char", "mode",
                    "median_time", "best_score", "median_replicates")],
      row.names = FALSE)

# Save to CSV
output_file <- if (length(commandArgs(TRUE)) > 0) {
  commandArgs(TRUE)[1]
} else {
  sprintf("tests/benchmark/bench-results-%s.csv", format(Sys.Date()))
}
write.csv(bench_df, output_file, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", output_file))
