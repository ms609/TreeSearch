#!/usr/bin/env Rscript
# Performance regression benchmark for TreeSearch C++ engine.
# Run after every significant code change to catch quality or speed regressions.
#
# Usage:
#   Rscript inst/benchmarks/bench_regression.R [lib_path]
#   Rscript inst/benchmarks/bench_regression.R --datasets=Vinther2008,Zhu2013 --budget=30
#   Rscript inst/benchmarks/bench_regression.R --datasets=all --budget=20 --output=results.csv
#
# Arguments (positional, legacy):
#   lib_path    Library path for TreeSearch (default: auto-detect)
#
# Arguments (named):
#   --lib=PATH        Library path (overrides positional)
#   --datasets=NAMES  Comma-separated dataset names, or "all" (default: core 3)
#   --budget=SECS     Per-dataset time budget in seconds (default: 30)
#   --output=FILE     Write CSV results to FILE (in addition to stdout)
#   --threads=N       Number of threads (default: 1)
#
# Each benchmark runs in its own subprocess to isolate any crashes.
#
# Asserts:
#   1. Score quality: each dataset must reach its max allowed score.
#   2. Timing: no dataset should take more than 3x its reference time.
#
# Exit code 0 = pass, 1 = regression detected.

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

named_args <- list()
positional_args <- character(0)
for (arg in args) {
  if (grepl("^--", arg)) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    named_args[[parts[1]]] <- if (length(parts) > 1) parts[2] else "true"
  } else {
    positional_args <- c(positional_args, arg)
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a

lib_path <- named_args[["lib"]] %||%
  (if (length(positional_args)) positional_args[1] else NULL)
budget <- as.numeric(named_args[["budget"]] %||% "30")
output_file <- named_args[["output"]]
n_threads <- as.integer(named_args[["threads"]] %||% "1")
dataset_arg <- named_args[["datasets"]]

# --- Reference data ---
# Max scores are ~1-2% above optimal to allow for stochastic variation.
# ref_time_s is the expected time at budget=30s with 1 thread.
all_benchmarks <- list(
  Vinther2008 = list(n_tip = 23, max_score = 80, ref_time_s = 1.0),
  Agnarsson2004 = list(n_tip = 62, max_score = 785, ref_time_s = 5.0),
  Zhu2013 = list(n_tip = 75, max_score = 662, ref_time_s = 8.0),
  Longrich2010 = list(n_tip = 20, max_score = 132, ref_time_s = 0.5),
  Sansom2010 = list(n_tip = 23, max_score = 190, ref_time_s = 0.8),
  DeAssis2011 = list(n_tip = 33, max_score = 66, ref_time_s = 1.0),
  Aria2015 = list(n_tip = 35, max_score = 145, ref_time_s = 1.5),
  Wortley2006 = list(n_tip = 37, max_score = 500, ref_time_s = 2.0),
  Griswold1999 = list(n_tip = 43, max_score = 415, ref_time_s = 3.0),
  Schulze2007 = list(n_tip = 52, max_score = 168, ref_time_s = 4.0),
  Eklund2004 = list(n_tip = 54, max_score = 450, ref_time_s = 4.0),
  Zanol2014 = list(n_tip = 74, max_score = 1345, ref_time_s = 7.0),
  Giles2015 = list(n_tip = 78, max_score = 725, ref_time_s = 7.0),
  Dikow2009 = list(n_tip = 88, max_score = 1625, ref_time_s = 10.0)
)

# Select datasets
default_names <- c("Vinther2008", "Agnarsson2004", "Zhu2013")
if (is.null(dataset_arg) || dataset_arg == "") {
  selected_names <- default_names
} else if (tolower(dataset_arg) == "all") {
  selected_names <- names(all_benchmarks)
} else {
  selected_names <- trimws(strsplit(dataset_arg, ",")[[1]])
  unknown <- setdiff(selected_names, names(all_benchmarks))
  if (length(unknown)) {
    stop("Unknown datasets: ", paste(unknown, collapse = ", "),
         "\nAvailable: ", paste(names(all_benchmarks), collapse = ", "))
  }
}

benchmarks <- all_benchmarks[selected_names]

# Resolve library path
if (is.null(lib_path)) {
  candidates <- c(Sys.glob(".agent-*"), Sys.glob(".builds/TreeSearch-*"))
  if (length(candidates)) {
    lib_path <- candidates[1]
    cat("Auto-detected library:", lib_path, "\n")
  } else {
    lib_path <- .libPaths()[1]
  }
}

cat("=== TreeSearch Performance Regression Benchmark ===\n")
cat(sprintf("  Library:  %s\n", lib_path))
cat(sprintf("  Datasets: %s\n", paste(selected_names, collapse = ", ")))
cat(sprintf("  Budget:   %ds per dataset\n", budget))
cat(sprintf("  Threads:  %d\n\n", n_threads))

n_pass <- 0L
n_fail <- 0L
results <- list()

for (nm in names(benchmarks)) {
  bm <- benchmarks[[nm]]
  cat(sprintf("--- %s (%d tips) ---\n", nm, bm$n_tip))

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
      maxReplicates = 100L, targetHits = 3L,
      ratchetCycles = 12L, driftCycles = 2L,
      xssRounds = 3L, xssPartitions = 4L,
      rssRounds = 1L, cssRounds = 0L,
      fuseInterval = 3L, maxSeconds = %d,
      verbosity = 0L, nThreads = %dL
    )
    elapsed <- (proc.time() - t0)[3]
    cat(result$best_score, elapsed, result$replicates, sep = " ")
  ', lib_path, nm, budget, n_threads)

  tf <- tempfile(fileext = ".R")
  writeLines(script, tf)
  timeout <- max(budget * 3, 60)
  output <- tryCatch(
    system2("Rscript", tf, stdout = TRUE, stderr = FALSE, timeout = timeout),
    error = function(e) paste("ERROR:", conditionMessage(e))
  )
  unlink(tf)

  if (length(output) == 0 || startsWith(output[length(output)], "ERROR")) {
    cat("  CRASHED or timed out\n")
    n_fail <- n_fail + 1L
    results[[nm]] <- data.frame(
      dataset = nm, n_tip = bm$n_tip,
      score = NA, elapsed = NA, replicates = NA, status = "CRASH",
      stringsAsFactors = FALSE
    )
    next
  }

  vals <- strsplit(trimws(output[length(output)]), "\s+")[[1]]
  if (length(vals) < 2) {
    cat("  Unexpected output:", output[length(output)], "\n")
    n_fail <- n_fail + 1L
    results[[nm]] <- data.frame(
      dataset = nm, n_tip = bm$n_tip,
      score = NA, elapsed = NA, replicates = NA, status = "ERROR",
      stringsAsFactors = FALSE
    )
    next
  }

  score <- as.numeric(vals[1])
  elapsed <- as.numeric(vals[2])
  reps <- if (length(vals) >= 3) as.integer(vals[3]) else NA_integer_

  score_ok <- score <= bm$max_score
  time_limit <- bm$ref_time_s * 3 * (budget / 30)
  time_ok <- elapsed <= time_limit

  status <- if (score_ok && time_ok) "PASS" else "FAIL"
  cat(sprintf("  Score: %.0f (max: %d) %s\n",
              score, bm$max_score,
              if (score_ok) "OK" else "REGRESSION"))
  cat(sprintf("  Time:  %.2fs (limit: %.1fs) %s\n",
              elapsed, time_limit,
              if (time_ok) "OK" else "REGRESSION"))
  if (!is.na(reps)) cat(sprintf("  Reps:  %d\n", reps))
  cat(sprintf("  Result: %s\n\n", status))

  if (status == "PASS") n_pass <- n_pass + 1L
  else n_fail <- n_fail + 1L

  results[[nm]] <- data.frame(
    dataset = nm, n_tip = bm$n_tip,
    score = score, elapsed = elapsed, replicates = reps, status = status,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("=== Summary: %d PASS, %d FAIL ===\n", n_pass, n_fail))

# Write CSV output if requested
if (!is.null(output_file)) {
  df <- do.call(rbind, results)
  df$budget_s <- budget
  df$threads <- n_threads
  df$timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, output_file, row.names = FALSE)
  cat(sprintf("Results written to %s\n", output_file))
}

if (n_fail > 0L) {
  cat("\nREGRESSIONS DETECTED.\n")
  quit(status = 1L)
} else {
  cat("\nAll benchmarks passed.\n")
  quit(status = 0L)
}
