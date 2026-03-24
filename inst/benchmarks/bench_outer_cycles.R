# bench_outer_cycles.R
#
# Compares thorough preset with outerCycles=1 vs outerCycles=2 across all 14
# standard benchmark datasets. Uses 3 seeds x 20s time budget per condition.
#
# Run from package root via:
#   Rscript inst/benchmarks/bench_outer_cycles.R
#
# Results saved to inst/benchmarks/results_outer_cycles.csv

.libPaths(c(".agent-X", .libPaths()))
library(TreeSearch)
library(TreeTools)

SRC <- getwd()
source(file.path(SRC, "inst/benchmarks/bench_datasets.R"))
source(file.path(SRC, "inst/benchmarks/bench_framework.R"))

BUDGET_S  <- 20
SEEDS     <- c(1031L, 2847L, 7193L)
OUT_FILE  <- file.path(SRC, "inst/benchmarks/results_outer_cycles.csv")

cat("TreeSearch version:", as.character(packageVersion("TreeSearch")), "\n")
cat(sprintf("Budget: %ds | Seeds: %d\n", BUDGET_S, length(SEEDS)))

# Build thorough strategy base (matches get_strategy("thorough") in bench_framework.R)
thorough_base <- list(
  wagnerStarts          = 3L,
  tbrMaxHits            = 3L,
  tabuSize              = 200L,
  ratchetCycles         = 20L,
  ratchetPerturbProb    = 0.25,
  ratchetPerturbMode    = 2L,
  ratchetPerturbMaxMoves = 5L,
  ratchetAdaptive       = TRUE,
  driftCycles           = 12L,
  driftAfdLimit         = 5L,
  driftRfdLimit         = 0.15,
  xssRounds             = 5L,
  xssPartitions         = 6L,
  rssRounds             = 3L,
  cssRounds             = 2L,
  cssPartitions         = 6L,
  sectorMinSize         = 6L,
  sectorMaxSize         = 80L,
  fuseInterval          = 2L,
  fuseAcceptEqual       = TRUE,
  nniFirst              = TRUE,
  sprFirst              = FALSE,
  consensusStableReps   = 3L
)

conditions <- list(
  thorough_1 = c(thorough_base, list(outerCycles = 1L)),
  thorough_2 = c(thorough_base, list(outerCycles = 2L))
)

datasets <- load_benchmark_datasets()
cat("Datasets loaded:", length(datasets), "\n\n")

total_runs <- length(BENCHMARK_NAMES) * length(conditions) * length(SEEDS)
cat(sprintf("Total runs: %d x %d conditions x %d seeds = %d\n\n",
            length(BENCHMARK_NAMES), length(conditions), length(SEEDS), total_runs))

rows <- list()
idx  <- 0L

for (ds_name in BENCHMARK_NAMES) {
  ds <- datasets[[ds_name]]
  if (is.null(ds)) { warning("Skipping ", ds_name); next }

  for (cond_name in names(conditions)) {
    strat <- conditions[[cond_name]]

    for (seed in SEEDS) {
      idx <- idx + 1L
      cat(sprintf("[%3d/%d] %-14s | %-12s | seed %d ... ",
                  idx, total_runs, ds_name, cond_name, seed))
      flush.console()

      t_start <- proc.time()
      set.seed(seed)
      result <- tryCatch(
        do.call(TreeSearch:::ts_driven_search,
                c(list(contrast    = ds$contrast,
                        tip_data    = ds$tip_data,
                        weight      = ds$weight,
                        levels      = ds$levels,
                        maxReplicates = 200L,
                        targetHits  = max(10L, ds$n_taxa %/% 5L),
                        maxSeconds  = as.double(BUDGET_S),
                        verbosity   = 0L),
                  strat)),
        error = function(e) {
          cat("ERROR:", conditionMessage(e), "\n"); NULL
        }
      )
      wall_s <- as.double((proc.time() - t_start)[3])

      if (is.null(result)) {
        rows[[idx]] <- data.frame(
          dataset = ds_name, condition = cond_name, seed = seed,
          n_taxa = ds$n_taxa, best_score = NA_real_,
          replicates = NA_integer_, hits_to_best = NA_integer_,
          wall_s = wall_s, stringsAsFactors = FALSE
        )
        next
      }

      cat(sprintf("score=%.0f  reps=%d  wall=%.1fs\n",
                  result$best_score, result$replicates, wall_s))

      rows[[idx]] <- data.frame(
        dataset       = ds_name,
        condition     = cond_name,
        seed          = seed,
        n_taxa        = ds$n_taxa,
        best_score    = result$best_score,
        replicates    = result$replicates,
        hits_to_best  = result$hits_to_best,
        wall_s        = wall_s,
        stringsAsFactors = FALSE
      )
    }
  }
}

results_df <- do.call(rbind, rows)
write.csv(results_df, OUT_FILE, row.names = FALSE)
cat("\nResults written to:", OUT_FILE, "\n")

# Quick summary
library(dplyr)
summary_tbl <- results_df |>
  filter(!is.na(best_score)) |>
  group_by(dataset, n_taxa, condition) |>
  summarise(median_score = median(best_score),
            median_reps  = median(replicates),
            .groups = "drop") |>
  tidyr::pivot_wider(names_from = condition,
                     values_from = c(median_score, median_reps)) |>
  mutate(delta = median_score_thorough_2 - median_score_thorough_1) |>
  arrange(n_taxa)

cat("\n===== outerCycles=2 vs outerCycles=1 (lower score = better) =====\n")
cat(sprintf("%-16s %5s  %8s  %8s  %6s  %5s  %5s\n",
            "Dataset", "Tips", "OC1_score", "OC2_score", "Delta",
            "OC1_reps", "OC2_reps"))
cat(strrep("-", 68), "\n")
for (i in seq_len(nrow(summary_tbl))) {
  r <- summary_tbl[i, ]
  cat(sprintf("%-16s %5d  %8.0f  %8.0f  %+6.1f  %5.0f  %5.0f\n",
              r$dataset, r$n_taxa,
              r$median_score_thorough_1, r$median_score_thorough_2,
              r$delta,
              r$median_reps_thorough_1, r$median_reps_thorough_2))
}
improved  <- sum(summary_tbl$delta < -0.5, na.rm = TRUE)
unchanged <- sum(abs(summary_tbl$delta) <= 0.5, na.rm = TRUE)
worse     <- sum(summary_tbl$delta > 0.5, na.rm = TRUE)
cat(strrep("-", 68), "\n")
cat(sprintf("Improved: %d  Unchanged: %d  Worse: %d\n",
            improved, unchanged, worse))
