#!/usr/bin/env Rscript
# T-290: Brazeau-track baseline benchmark runner
#
# Runs the fixed 20-matrix Brazeau sample under EW + IW(k=10) with
# default and thorough strategy presets. Includes sample-size stability
# analysis on the EW results.
#
# Usage:
#   Rscript t290_run.R <output_dir>
#
# Requires: TreeSearch (installed), neotrans corpus in ../neotrans/

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."

repo_root <- getwd()
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) getwd()
  )
  repo_root <- normalizePath(file.path(script_dir, "..", ".."),
                             mustWork = FALSE)
}
setwd(repo_root)

cat("=== T-290: Brazeau-Track Baseline Benchmark ===\n")
cat("Repo root:", repo_root, "\n")
cat("Output dir:", outdir, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

library(TreeSearch)
library(TreeTools)

source("dev/benchmarks/bench_datasets.R")
source("dev/benchmarks/bench_brazeau_baseline.R")

# ---- Configuration ----
TIMEOUT <- 30
N_SEEDS <- 5L
STRATEGIES <- c("default", "thorough")

# ---- Run full grid ----
cat("Phase 1: Full benchmark grid\n")
cat(sprintf("  Datasets: %d (MBANK_BRAZEAU_SAMPLE)\n",
            length(MBANK_BRAZEAU_SAMPLE)))
cat(sprintf("  Weightings: EW + IW(k=10)\n"))
cat(sprintf("  Strategies: %s\n", paste(STRATEGIES, collapse = ", ")))
cat(sprintf("  Timeout: %ds, Seeds: %d\n\n", TIMEOUT, N_SEEDS))

results <- run_brazeau_grid(
  dataset_keys = MBANK_BRAZEAU_SAMPLE,
  weightings = list(EW = Inf, IW10 = 10),
  strategies = STRATEGIES,
  timeout_s = TIMEOUT,
  seeds = seq_len(N_SEEDS),
  hits = 10L,
  reps = 50L
)

# Save raw results
outfile <- file.path(outdir,
                     sprintf("t290_brazeau_%s.csv",
                             format(Sys.time(), "%Y%m%d_%H%M")))
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nRaw results saved: %s (%d rows)\n\n", outfile, nrow(results)))

# ---- Summary ----
cat("=== Summary ===\n")
summary_df <- summarize_brazeau(results)
print(summary_df)

outfile_summary <- file.path(outdir,
                             sprintf("t290_summary_%s.csv",
                                     format(Sys.time(), "%Y%m%d_%H%M")))
write.csv(summary_df, outfile_summary, row.names = FALSE)
cat(sprintf("Summary saved: %s\n\n", outfile_summary))

# ---- Phase 2: Sample-size stability (EW only) ----
cat("=== Phase 2: Sample-Size Stability Analysis (EW) ===\n")
ew_results <- results[results$weighting == "EW", ]

if (length(unique(ew_results$strategy)) >= 2) {
  stability <- sample_size_stability(ew_results, n_boot = 500L)
  cat("Subsample size vs ranking stability (Kendall's W):\n")
  print(stability)
  
  outfile_stability <- file.path(outdir,
                                 sprintf("t290_stability_%s.csv",
                                         format(Sys.time(), "%Y%m%d_%H%M")))
  write.csv(stability, outfile_stability, row.names = FALSE)
  cat(sprintf("Stability results saved: %s\n", outfile_stability))
  
  # Check if 20 is enough
  if (nrow(stability) >= 2) {
    last_two <- tail(stability, 2)
    delta_W <- diff(last_two$median_W)
    cat(sprintf("\nΔW at last two increments: %.3f\n", delta_W))
    if (abs(delta_W) < 0.02) {
      cat("Ranking has plateaued — 20 matrices appear sufficient.\n")
    } else {
      cat("Ranking still changing — consider expanding the sample.\n")
    }
  }
} else {
  cat("Skipped: need >= 2 strategies for stability analysis.\n")
}

cat(sprintf("\n=== Completed: %s ===\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
