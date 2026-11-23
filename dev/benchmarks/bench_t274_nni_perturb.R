# bench_t274_nni_perturb.R
#
# T-274: Benchmark nniPerturbCycles=0 vs 5 at thorough-preset scale.
#
# S-PROF round 6 found NNI-perturb = 34.3% of Zhu2013 (75t) thorough-preset
# search time with only 14% hit rate and ~1-step mean improvement.
# This benchmark tests whether removing NNI-perturb improves time-adjusted
# expected best score at 30s and 60s budgets on 65–88 tip datasets.
#
# METHODOLOGY: Per-replicate sampling.
#   - maxReplicates=1 per run, many seeds → per-replicate score distribution
#   - time_per_rep estimated from wall time
#   - expected_best(scores, k=floor(budget/median_time)) at 30s/60s
#
# Usage:
#   Rscript dev/benchmarks/bench_t274_nni_perturb.R [lib_path]
#   Default lib_path = .agent-F
#
# Results: dev/benchmarks/results_t274_nni_perturb.csv
# Run time: ~12-18 min (3 datasets x 2 conditions x 20 seeds)

args <- commandArgs(trailingOnly = TRUE)
lib_path <- if (length(args) >= 1) args[[1L]] else ".agent-F"
.libPaths(c(lib_path, .libPaths()))
library(TreeSearch)
library(TreeTools)

cat("TreeSearch version:", as.character(packageVersion("TreeSearch")), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
DATASETS  <- c("Zhu2013", "Giles2015", "Dikow2009")  # 75, 78, 88 tips
BUDGETS_S <- c(30, 60)
N_SEEDS   <- 20L
NNI_CONDITIONS <- c(0L, 5L)
OUT_FILE  <- "dev/benchmarks/results_t274_nni_perturb.csv"

# Seeds — fixed for reproducibility
set.seed(4718)
seeds <- sample.int(99999L, N_SEEDS)

# ------------------------------------------------------------
# expected_best: bootstrap estimate of expected minimum from k draws
# ------------------------------------------------------------
expected_best <- function(scores, k, n_boot = 5000L) {
  mean(replicate(n_boot, min(sample(scores, k, replace = TRUE))))
}

# ------------------------------------------------------------
# Per-replicate runs
# ------------------------------------------------------------
total_runs <- length(DATASETS) * length(NNI_CONDITIONS) * N_SEEDS
cat(sprintf("Total runs: %d datasets x %d conditions x %d seeds = %d\n\n",
            length(DATASETS), length(NNI_CONDITIONS), N_SEEDS, total_runs))

rows <- list()
idx  <- 0L

for (ds_name in DATASETS) {
  dataset <- TreeSearch::inapplicable.phyData[[ds_name]]
  if (is.null(dataset)) {
    warning("Dataset not found: ", ds_name)
    next
  }
  n_taxa <- length(dataset)
  n_char <- sum(attr(dataset, "weight"))

  cat(sprintf("=== %s (%dt, %dc) ===\n", ds_name, n_taxa, n_char))

  for (nni_cycles in NNI_CONDITIONS) {
    cond_label <- if (nni_cycles == 0L) "nni=0" else sprintf("nni=%d", nni_cycles)

    for (seed in seeds) {
      idx <- idx + 1L
      cat(sprintf("[%3d/%d] %-12s | %-6s | seed %5d ... ",
                  idx, total_runs, ds_name, cond_label, seed))
      flush.console()

      set.seed(seed)
      t0 <- proc.time()[[3L]]
      result <- tryCatch(
        # Pass nniPerturbCycles via ... so it overrides the thorough preset
        # for just that parameter, leaving all other thorough params intact.
        MaximizeParsimony(dataset,
                          strategy       = "thorough",
                          nniPerturbCycles = as.integer(nni_cycles),
                          maxReplicates  = 1L,
                          nThreads       = 1L,
                          verbosity      = 0L),
        error = function(e) {
          cat("ERROR:", conditionMessage(e), "\n")
          NULL
        }
      )
      wall_s <- proc.time()[[3L]] - t0

      if (is.null(result)) {
        rows[[idx]] <- data.frame(
          dataset = ds_name, n_taxa = n_taxa, nni_cycles = nni_cycles,
          seed = seed, best_score = NA_real_, wall_s = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      best_score <- min(attr(result, "score"), na.rm = TRUE)
      cat(sprintf("score=%.0f  wall=%.1fs\n", best_score, wall_s))

      rows[[idx]] <- data.frame(
        dataset    = ds_name,
        n_taxa     = n_taxa,
        nni_cycles = nni_cycles,
        seed       = seed,
        best_score = best_score,
        wall_s     = wall_s,
        stringsAsFactors = FALSE
      )
    }
    cat("\n")
  }
}

results_df <- do.call(rbind, rows)
write.csv(results_df, OUT_FILE, row.names = FALSE)
cat("\nResults written to:", OUT_FILE, "\n\n")

# ------------------------------------------------------------
# Analysis: Time-adjusted expected best
# ------------------------------------------------------------
cat("===== Time-adjusted expected best (lower score = better) =====\n\n")

for (ds_name in DATASETS) {
  sub <- results_df[results_df$dataset == ds_name & !is.na(results_df$best_score), ]
  cat(sprintf("--- %s ---\n", ds_name))

  for (budget in BUDGETS_S) {
    cat(sprintf("  Budget = %ds:\n", budget))
    for (nni in NNI_CONDITIONS) {
      d <- sub[sub$nni_cycles == nni, ]
      if (nrow(d) < 5L) { cat(sprintf("    nni=%d: insufficient data\n", nni)); next }
      med_time <- median(d$wall_s, na.rm = TRUE)
      k <- max(1L, floor(budget / med_time))
      eb <- expected_best(d$best_score, k)
      cat(sprintf("    nni=%d: median_time=%.1fs, k=%d reps, expected_best=%.1f  (n=%d)\n",
                  nni, med_time, k, eb, nrow(d)))
    }
  }
  cat("\n")
}

# Summary table: delta (nni=0 - nni=5) at each budget
cat("===== Expected-best delta (nni=0 vs nni=5, negative = nni=0 better) =====\n")
cat(sprintf("%-14s  %8s  %8s  %8s  %8s\n",
            "Dataset", "30s_nni0", "30s_nni5", "60s_nni0", "60s_nni5"))
cat(strrep("-", 56), "\n")

for (ds_name in DATASETS) {
  sub <- results_df[results_df$dataset == ds_name & !is.na(results_df$best_score), ]
  row_vals <- c(ds_name)

  for (budget in BUDGETS_S) {
    for (nni in NNI_CONDITIONS) {
      d <- sub[sub$nni_cycles == nni, ]
      if (nrow(d) < 5L) { row_vals <- c(row_vals, "N/A"); next }
      med_time <- median(d$wall_s, na.rm = TRUE)
      k <- max(1L, floor(budget / med_time))
      eb <- expected_best(d$best_score, k)
      row_vals <- c(row_vals, sprintf("%.1f", eb))
    }
  }
  cat(sprintf("%-14s  %8s  %8s  %8s  %8s\n",
              row_vals[[1L]], row_vals[[2L]], row_vals[[3L]],
              row_vals[[4L]], row_vals[[5L]]))
}

cat("\n")
cat("Interpretation: Positive delta = nni=0 is better (removes overhead).\n")
cat("Negative delta = nni=5 is better (perturbation value exceeds overhead).\n")
