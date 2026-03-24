# bench_large_preset.R
#
# Validates the T-179 "large" strategy preset against "thorough" on the
# 180-taxon mbank_X30754 dataset.
#
# Run from package root:
#   Rscript inst/benchmarks/bench_large_preset.R
#
# Results saved to inst/benchmarks/results_large_preset.csv

.libPaths(c(".agent-X", .libPaths()))
library(TreeSearch)
library(TreeTools)

SRC <- getwd()
source(file.path(SRC, "inst/benchmarks/bench_datasets.R"))
# Pull updated presets from source (no rebuild needed for pure-R changes)
source(file.path(SRC, "R/SearchControl.R"))
source(file.path(SRC, "R/MaximizeParsimony.R"))

BUDGET_S <- 60          # 60s per run — allows ~1 replicate at 180 tips
SEEDS    <- c(1031L, 2847L, 7193L, 4561L, 8822L)
OUT_FILE <- file.path(SRC, "inst/benchmarks/results_large_preset.csv")

cat("TreeSearch version:", as.character(packageVersion("TreeSearch")), "\n")
cat(sprintf("Budget: %ds | Seeds: %d\n\n", BUDGET_S, length(SEEDS)))

# Load 180-taxon dataset
large_ds_list <- load_large_benchmark_datasets()
ds_180 <- large_ds_list[["mbank_X30754"]]
if (is.null(ds_180)) stop("mbank_X30754 not found")
cat(sprintf("Dataset: mbank_X30754 | %d taxa | %d patterns\n\n",
            ds_180$n_taxa, length(ds_180$weight)))

# Use R-level SearchControl presets (sourced above)
presets <- .StrategyPresets()
conditions <- list(
  large    = unclass(presets[["large"]]),
  thorough = unclass(presets[["thorough"]])
)
conditions <- lapply(conditions, function(x) { attr(x, "class") <- NULL; x })

total_runs <- length(conditions) * length(SEEDS)
cat(sprintf("Total runs: %d conditions x %d seeds = %d\n\n",
            length(conditions), length(SEEDS), total_runs))

rows <- list()
idx  <- 0L

for (cond_name in names(conditions)) {
  strat <- conditions[[cond_name]]
  for (seed in SEEDS) {
    idx <- idx + 1L
    cat(sprintf("[%d/%d] %-10s | seed %d ... ",
                idx, total_runs, cond_name, seed))
    flush.console()

    t_start <- proc.time()
    set.seed(seed)
    result <- tryCatch(
      do.call(TreeSearch:::ts_driven_search,
              c(list(contrast    = ds_180$contrast,
                      tip_data    = ds_180$tip_data,
                      weight      = ds_180$weight,
                      levels      = ds_180$levels,
                      maxReplicates = 500L,
                      targetHits  = max(10L, ds_180$n_taxa %/% 5L),
                      maxSeconds  = as.double(BUDGET_S),
                      verbosity   = 0L),
                strat)),
      error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
    )
    wall_s <- as.double((proc.time() - t_start)[3])

    if (is.null(result)) next

    cat(sprintf("score=%.0f  reps=%d  wall=%.1fs\n",
                result$best_score, result$replicates, wall_s))

    rows[[idx]] <- data.frame(
      condition = cond_name, seed = seed,
      best_score = result$best_score,
      replicates = result$replicates,
      hits_to_best = result$hits_to_best,
      wall_s = wall_s,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, rows)
write.csv(results_df, OUT_FILE, row.names = FALSE)
cat("\nResults written to:", OUT_FILE, "\n")

# Summary
cat("\n===== large vs thorough on mbank_X30754 (180 tips, 60s budget) =====\n")
cat(sprintf("%-12s  %8s  %8s  %8s  %8s\n",
            "Condition", "Min", "Median", "Max", "Med.reps"))
for (cond in names(conditions)) {
  r <- results_df[results_df$condition == cond & !is.na(results_df$best_score), ]
  cat(sprintf("%-12s  %8.0f  %8.0f  %8.0f  %8.0f\n",
              cond, min(r$best_score), median(r$best_score),
              max(r$best_score), median(r$replicates)))
}

# Per-seed comparison
cat("\nPer-seed comparison (large - thorough, negative = large better):\n")
for (s in SEEDS) {
  lrg <- results_df$best_score[results_df$condition == "large" & results_df$seed == s]
  thr <- results_df$best_score[results_df$condition == "thorough" & results_df$seed == s]
  if (length(lrg) == 1 && length(thr) == 1) {
    cat(sprintf("  seed %d: large=%4.0f  thorough=%4.0f  delta=%+.0f\n",
                s, lrg, thr, lrg - thr))
  }
}
