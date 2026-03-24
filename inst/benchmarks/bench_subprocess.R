# Subprocess-isolated benchmark: each run in its own Rscript process.
# Workaround for T-025 (ratchet-triggered optimization-dependent UB that
# causes segfaults on consecutive ts_driven_search calls).

library(TreeSearch)
library(TreeTools)

source("inst/benchmarks/bench_datasets.R")
source("inst/benchmarks/bench_framework.R")

GRID_DATASETS <- c(
  "Longrich2010",   # 20 tips
  "Vinther2008",    # 23 tips
  "Aria2015",       # 35 tips
  "Griswold1999",   # 43 tips
  "Agnarsson2004",  # 62 tips
  "Zhu2013",        # 75 tips
  "Giles2015",      # 78 tips
  "Dikow2009"       # 88 tips
)

# One benchmark in a subprocess; returns CSV line or NA on crash
run_one_subprocess <- function(ds_name, strat_name, seed, maxSeconds = 20L,
                                maxReplicates = 100L) {
  script <- sprintf('
library(TreeSearch, lib.loc = if (dir.exists(".agent-a")) ".agent-a" else .libPaths())
library(TreeTools)
source("inst/benchmarks/bench_datasets.R")
source("inst/benchmarks/bench_framework.R")
ds <- prepare_ts_data(TreeSearch::inapplicable.phyData[["%s"]])
strat <- get_strategy("%s")
targetHits <- max(10L, ds$n_taxa %%/%%  5L)
args <- c(
  list(contrast = ds$contrast, tip_data = ds$tip_data,
       weight = ds$weight, levels = ds$levels,
       maxReplicates = %dL, targetHits = targetHits,
       maxSeconds = %d, verbosity = 0L),
  strat)
set.seed(%dL)
t0 <- proc.time()
result <- do.call(TreeSearch:::ts_driven_search, args)
wall <- as.double((proc.time() - t0)[3])
cat(result$best_score, result$replicates, result$hits_to_best,
    result$pool_size, as.integer(result$timed_out), wall,
    result$timings[["wagner_ms"]], result$timings[["tbr_ms"]],
    result$timings[["xss_ms"]], result$timings[["rss_ms"]],
    result$timings[["css_ms"]], result$timings[["ratchet_ms"]],
    result$timings[["drift_ms"]], result$timings[["final_tbr_ms"]],
    result$timings[["fuse_ms"]], sep = ",")
', ds_name, strat_name, maxReplicates, maxSeconds, seed)

  tf <- tempfile(fileext = ".R")
  writeLines(script, tf)
  on.exit(unlink(tf))

  out <- tryCatch(
    system2("Rscript", c("--no-save", tf),
            stdout = TRUE, stderr = FALSE, timeout = maxSeconds + 30L),
    error = function(e) NA_character_
  )

  if (length(out) == 0 || is.na(out[1])) return(NULL)
  vals <- as.numeric(strsplit(out[length(out)], ",")[[1]])
  if (length(vals) != 15) return(NULL)

  data.frame(
    dataset = ds_name, strategy = strat_name, seed = seed,
    n_taxa = length(TreeSearch::inapplicable.phyData[[ds_name]]),
    best_score = vals[1], replicates = vals[2], hits_to_best = vals[3],
    pool_size = vals[4], timed_out = as.logical(vals[5]),
    wall_s = vals[6],
    wagner_ms = vals[7], tbr_ms = vals[8], xss_ms = vals[9],
    rss_ms = vals[10], css_ms = vals[11], ratchet_ms = vals[12],
    drift_ms = vals[13], final_tbr_ms = vals[14], fuse_ms = vals[15],
    stringsAsFactors = FALSE
  )
}

# Run grid using subprocess isolation
run_grid_safe <- function(dataset_names = GRID_DATASETS,
                           strategy_names = STRATEGY_NAMES,
                           replicates = 3L,
                           maxSeconds = 20L,
                           base_seed = 7142L) {
  n_combos <- length(dataset_names) * length(strategy_names) * replicates
  cat(sprintf("Grid: %d datasets x %d strategies x %d reps = %d runs (subprocess)\n",
              length(dataset_names), length(strategy_names), replicates, n_combos))

  rows <- vector("list", n_combos)
  idx <- 0L

  for (ds_name in dataset_names) {
    for (strat_name in strategy_names) {
      for (rep in seq_len(replicates)) {
        idx <- idx + 1L
        seed <- base_seed + (idx - 1L) * 7L

        cat(sprintf("[%3d/%d] %-15s x %-16s rep %d ... ",
                    idx, n_combos, ds_name, strat_name, rep))

        res <- run_one_subprocess(ds_name, strat_name, seed,
                                   maxSeconds = maxSeconds)
        if (is.null(res)) {
          cat("CRASH/ERROR\n")
          next
        }
        cat(sprintf("score=%.0f wall=%.1fs reps=%d %s\n",
                    res$best_score, res$wall_s, res$replicates,
                    if (res$timed_out) "[TIMEOUT]" else ""))
        rows[[idx]] <- res
      }
    }
  }

  result <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  outfile <- "inst/benchmarks/results_grid.csv"
  write.csv(result, outfile, row.names = FALSE)
  cat(sprintf("\nResults saved to %s (%d rows)\n", outfile, nrow(result)))
  invisible(result)
}

# Main
cat("Starting subprocess-isolated benchmark grid...\n\n")
results <- run_grid_safe()
