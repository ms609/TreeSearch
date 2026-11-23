#!/usr/bin/env Rscript
# T-252: MorphoBank training-set baseline benchmark
#
# Runs the fixed 25-matrix training sample at 30s, 60s, and 120s budgets
# with the "default" strategy, 5 seeds per combination.
# Total: 25 matrices x 3 budgets x 5 seeds = 375 runs.
# Estimated wall time: ~4–5 hours (most runs hit timeout).
#
# Usage:
#   Rscript bench_t252_mbank_training.R <output_dir>
#
# Requires: TreeSearch (installed), neotrans corpus in ../neotrans/

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."

# Find the repo root (this script lives in dev/benchmarks/)
# When run from repo root (cd $REPO; Rscript dev/benchmarks/...), getwd() is it.
repo_root <- getwd()
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  # Try relative to script location
  script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) getwd()
  )
  repo_root <- normalizePath(file.path(script_dir, "..", ".."),
                             mustWork = FALSE)
}
setwd(repo_root)

cat("=== T-252: MorphoBank Training-Set Benchmark ===\n")
cat("Repo root:", repo_root, "\n")
cat("Output dir:", outdir, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

library(TreeSearch)
library(TreeTools)

source("dev/benchmarks/bench_datasets.R")
source("dev/benchmarks/bench_framework.R")

# ---- Configuration ----
BUDGETS <- c(30, 60, 120)  # seconds
N_SEEDS <- 5L
BASE_SEED <- 3847L
STRATEGY <- "default"

# ---- Load training matrices ----
cat("Loading MorphoBank catalogue...\n")
catalogue <- load_mbank_catalogue()
cat(sprintf("Catalogue: %d usable matrices\n", nrow(catalogue)))

cat(sprintf("Loading %d fixed training matrices...\n",
            length(MBANK_FIXED_SAMPLE)))
datasets <- load_mbank_datasets(catalogue, keys = MBANK_FIXED_SAMPLE)
cat(sprintf("Successfully loaded: %d matrices\n\n", length(datasets)))

if (length(datasets) == 0) {
  stop("No datasets loaded. Is the neotrans repo available?")
}

# ---- Characterize datasets ----
cat("Dataset characteristics:\n")
chars <- do.call(rbind, lapply(names(datasets), function(nm) {
  ch <- characterize_dataset(datasets[[nm]])
  ch$key <- nm
  ch
}))
chars <- chars[order(chars$n_taxa), ]
print(chars[, c("key", "n_taxa", "n_chars", "n_patterns",
                "pct_missing", "pct_inapp", "n_app_states")])
cat("\n")

# ---- Run benchmarks ----
strat <- get_strategy(STRATEGY)
all_results <- list()

for (budget in BUDGETS) {
  cat(sprintf("\n========== Budget: %ds ==========\n", budget))

  results <- run_benchmark_grid(
    dataset_names = names(datasets),
    strategy_names = STRATEGY,
    replicates = N_SEEDS,
    maxReplicates = 100L,
    maxSeconds = budget,
    base_seed = BASE_SEED,
    datasets = datasets
  )
  results$budget_s <- budget
  results$source <- "mbank_training"

  # Save intermediate results per budget
  outfile <- file.path(
    outdir,
    sprintf("t252_mbank_%ds_%s.csv", budget,
            format(Sys.time(), "%Y%m%d_%H%M"))
  )
  write.csv(results, outfile, row.names = FALSE)
  cat(sprintf("Saved %d rows to %s\n", nrow(results), outfile))

  all_results[[as.character(budget)]] <- results
}

# ---- Combine and save final results ----
final <- do.call(rbind, all_results)
final_file <- file.path(outdir,
                        sprintf("t252_mbank_all_%s.csv",
                                format(Sys.time(), "%Y%m%d_%H%M")))
write.csv(final, final_file, row.names = FALSE)
cat(sprintf("\n=== Final results: %d rows saved to %s ===\n",
            nrow(final), final_file))

# ---- Summary statistics ----
cat("\n=== Summary by budget ===\n")
for (budget in BUDGETS) {
  sub <- final[final$budget_s == budget, ]
  cat(sprintf("\n--- %ds budget (%d runs) ---\n", budget, nrow(sub)))
  cat(sprintf("  Median score: %.1f\n", median(sub$best_score, na.rm = TRUE)))
  cat(sprintf("  Timed out: %d/%d (%.0f%%)\n",
              sum(sub$timed_out, na.rm = TRUE), nrow(sub),
              100 * mean(sub$timed_out, na.rm = TRUE)))
  cat(sprintf("  Median replicates: %.0f\n",
              median(sub$replicates, na.rm = TRUE)))
  cat(sprintf("  Median wall time: %.1fs\n",
              median(sub$wall_s, na.rm = TRUE)))
}

cat("\n=== Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
