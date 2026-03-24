#!/usr/bin/env Rscript
# Performance regression benchmark for TreeSearch C++ engine.
# Run after every significant code change to catch quality or speed regressions.
#
# Usage: Rscript inst/benchmarks/bench_regression.R [lib_path]
#
# Each benchmark runs in its own subprocess to isolate any crashes.
#
# Asserts:
#   1. Score quality: each dataset must reach its max allowed score.
#   2. Timing: no dataset should take more than 3x its reference time.
#
# Exit code 0 = pass, 1 = regression detected.

args <- commandArgs(trailingOnly = TRUE)
lib_path <- if (length(args) >= 1) args[1] else ".agent-c"

# --- Reference data ---
# Max scores are ~1-2% above optimal to allow for stochastic variation.
benchmarks <- list(
  list(
    name = "Vinther2008",
    n_tip = 23,
    max_score = 80,   # optimal 79
    ref_time_s = 1.0
  ),
  list(
    name = "Agnarsson2004",
    n_tip = 62,
    max_score = 785,  # optimal 778
    ref_time_s = 5.0
  ),
  list(
    name = "Zhu2013",
    n_tip = 75,
    max_score = 662,  # optimal 656
    ref_time_s = 8.0
  )
)

cat("=== TreeSearch Performance Regression Benchmark ===\n\n")

n_pass <- 0L
n_fail <- 0L

for (bm in benchmarks) {
  cat(sprintf("--- %s (%d tips) ---\n", bm$name, bm$n_tip))

  # Build subprocess script
  script <- sprintf('
    library(TreeSearch, lib.loc = "%s")
    library(TreeTools)
    ds <- TreeSearch::inapplicable.phyData[["%s"]]
    at <- attributes(ds)
    contrast <- at$contrast
    tip_data <- matrix(unlist(ds, use.names = FALSE), nrow = length(ds), byrow = TRUE)
    weight <- at$weight
    levels <- at$levels
    set.seed(4217)
    t0 <- proc.time()
    result <- TreeSearch:::ts_driven_search(
      contrast, tip_data, weight, levels,
      maxReplicates = 3L, targetHits = 2L,
      ratchetCycles = 5L, driftCycles = 2L,
      xssRounds = 3L, xssPartitions = 4L,
      rssRounds = 1L, cssRounds = 0L,
      fuseInterval = 3L, maxSeconds = 30,
      verbosity = 0L, nThreads = 1L
    )
    elapsed <- (proc.time() - t0)[3]
    cat(result$best_score, elapsed, sep = " ")
  ', lib_path, bm$name)

  tf <- tempfile(fileext = ".R")
  writeLines(script, tf)
  output <- tryCatch(
    system2("Rscript", tf, stdout = TRUE, stderr = FALSE, timeout = 60),
    error = function(e) paste("ERROR:", conditionMessage(e))
  )
  unlink(tf)

  if (length(output) == 0 || startsWith(output[length(output)], "ERROR")) {
    cat("  CRASHED or timed out\n")
    n_fail <- n_fail + 1L
    next
  }

  vals <- strsplit(trimws(output[length(output)]), "\\s+")[[1]]
  if (length(vals) < 2) {
    cat("  Unexpected output:", output[length(output)], "\n")
    n_fail <- n_fail + 1L
    next
  }

  score <- as.numeric(vals[1])
  elapsed <- as.numeric(vals[2])

  score_ok <- score <= bm$max_score
  time_ok <- elapsed <= bm$ref_time_s * 3

  status <- if (score_ok && time_ok) "PASS" else "FAIL"
  cat(sprintf("  Score: %.0f (max: %d) %s\n",
              score, bm$max_score,
              if (score_ok) "OK" else "REGRESSION"))
  cat(sprintf("  Time:  %.2fs (ref: %.1fs, limit: %.1fs) %s\n",
              elapsed, bm$ref_time_s, bm$ref_time_s * 3,
              if (time_ok) "OK" else "REGRESSION"))
  cat(sprintf("  Result: %s\n\n", status))

  if (status == "PASS") n_pass <- n_pass + 1L
  else n_fail <- n_fail + 1L
}

cat(sprintf("=== Summary: %d PASS, %d FAIL ===\n", n_pass, n_fail))

if (n_fail > 0L) {
  cat("\nREGRESSIONS DETECTED.\n")
  quit(status = 1L)
} else {
  cat("\nAll benchmarks passed.\n")
  quit(status = 0L)
}
