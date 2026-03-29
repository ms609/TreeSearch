# Experiment: TBR <-> XSS cycling under EW vs IW
#
# Question: After TBR converges, does XSS find improvements more often
# under implied weights than under equal weights?
#
# Protocol per (dataset x scoring x seed):
#   1. Random Wagner -> NNI -> TBR convergence
#   2. XSS on TBR-converged tree
#   3. TBR on XSS-modified tree
#   4. Repeat until neither improves (joint convergence), max 10 cycles
#
# Usage:
#   Rscript dev/benchmarks/expt_tbr_xss_cycling.R
#   # or source() interactively

library(TreeSearch)
library(TreeTools)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight   = at$weight,
    levels   = at$levels
  )
}

min_steps_for <- function(ds) {
  as.integer(TreeSearch::MinimumLength(ds, compress = TRUE))
}

concavity_value <- function(mode) {
  switch(mode,
    EW   = Inf,
    IW10 = 10,
    IW3  = 3,
    stop("Unknown scoring mode: ", mode)
  )
}

# Run a single Wagner -> NNI -> TBR -> [XSS <-> TBR cycling] experiment
run_one <- function(dataset, ds_data, min_steps, scoring_mode, seed,
                    max_cycles = 10L) {
  set.seed(seed)
  k <- concavity_value(scoring_mode)

  # Wagner start
  wagner <- TreeSearch:::ts_random_wagner_tree(
    ds_data$contrast, ds_data$tip_data, ds_data$weight, ds_data$levels,
    min_steps = min_steps, concavity = k
  )

  # NNI warmup
  nni <- TreeSearch:::ts_nni_search(
    wagner$edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
    ds_data$levels, maxHits = 20L, min_steps = min_steps, concavity = k
  )

  # TBR convergence (cycle 0)
  tbr <- TreeSearch:::ts_tbr_search(
    nni$edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
    ds_data$levels, maxHits = 1L, min_steps = min_steps, concavity = k
  )

  results <- data.frame(
    cycle            = 0L,
    phase            = "TBR",
    score_before     = nni$score,
    score_after      = tbr$score,
    improved         = tbr$score < nni$score - 1e-9,
    n_sectors_searched  = NA_integer_,
    n_sectors_improved  = NA_integer_,
    stringsAsFactors = FALSE
  )

  current_edge  <- tbr$edge
  current_score <- tbr$score

  for (cyc in seq_len(max_cycles)) {
    # XSS phase
    xss <- TreeSearch:::ts_xss_search(
      current_edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
      ds_data$levels,
      nPartitions = 4L, xssRounds = 3L, acceptEqual = FALSE,
      maxHits = 1L, min_steps = min_steps, concavity = k
    )

    xss_improved <- xss$score < current_score - 1e-9

    results <- rbind(results, data.frame(
      cycle            = cyc,
      phase            = "XSS",
      score_before     = current_score,
      score_after      = xss$score,
      improved         = xss_improved,
      n_sectors_searched  = xss$n_sectors_searched,
      n_sectors_improved  = xss$n_sectors_improved,
      stringsAsFactors = FALSE
    ))

    current_edge  <- xss$edge
    current_score <- xss$score

    # TBR phase
    tbr2 <- TreeSearch:::ts_tbr_search(
      current_edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
      ds_data$levels, maxHits = 1L, min_steps = min_steps, concavity = k
    )

    tbr_improved <- tbr2$score < current_score - 1e-9

    results <- rbind(results, data.frame(
      cycle            = cyc,
      phase            = "TBR",
      score_before     = current_score,
      score_after      = tbr2$score,
      improved         = tbr_improved,
      n_sectors_searched  = NA_integer_,
      n_sectors_improved  = NA_integer_,
      stringsAsFactors = FALSE
    ))

    current_edge  <- tbr2$edge
    current_score <- tbr2$score

    # Stop if neither improved
    if (!xss_improved && !tbr_improved) break
  }

  results$scoring <- scoring_mode
  results$seed    <- seed
  results
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

run_experiment <- function(n_seeds = 20L) {

  dataset_names <- c("Agnarsson2004", "Zhu2013", "Giles2015", "Dikow2009")
  scoring_modes <- c("EW", "IW10", "IW3")
  seeds <- sample.int(10000, n_seeds)

  all_results <- vector("list",
    length(dataset_names) * length(scoring_modes) * n_seeds)
  idx <- 0L

  for (ds_name in dataset_names) {
    message("=== ", ds_name, " ===")
    dataset   <- TreeSearch::inapplicable.phyData[[ds_name]]
    ds_data   <- make_ts_data(dataset)
    min_steps <- min_steps_for(dataset)

    for (mode in scoring_modes) {
      message("  ", mode, ": ", appendLF = FALSE)
      for (s in seq_along(seeds)) {
        if (s %% 5 == 0) message(s, " ", appendLF = FALSE)
        idx <- idx + 1L
        res <- run_one(dataset, ds_data, min_steps, mode, seeds[s])
        res$dataset <- ds_name
        all_results[[idx]] <- res
      }
      message()
    }
  }

  do.call(rbind, all_results)
}

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

analyse_results <- function(df) {

  # Focus on XSS steps after the initial TBR convergence (cycle >= 1)
  xss <- df[df$phase == "XSS" & df$cycle >= 1, ]

  message("\n--- XSS improvement rate after TBR convergence ---")
  agg <- aggregate(improved ~ dataset + scoring, data = xss, FUN = mean)
  names(agg)[3] <- "improvement_rate"
  print(agg[order(agg$dataset, agg$scoring), ], row.names = FALSE)

  message("\n--- Mean XSS improvement magnitude (when improved) ---")
  xss_improved <- xss[xss$improved, ]
  if (nrow(xss_improved) > 0) {
    xss_improved$delta <- xss_improved$score_before - xss_improved$score_after
    agg2 <- aggregate(delta ~ dataset + scoring, data = xss_improved, FUN = mean)
    names(agg2)[3] <- "mean_improvement"
    print(agg2[order(agg2$dataset, agg2$scoring), ], row.names = FALSE)
  } else {
    message("  (no XSS improvements found)")
  }

  message("\n--- Cycles to joint TBR+XSS convergence ---")
  conv <- aggregate(cycle ~ dataset + scoring + seed, data = df, FUN = max)
  agg3 <- aggregate(cycle ~ dataset + scoring, data = conv, FUN = mean)
  names(agg3)[3] <- "mean_cycles"
  print(agg3[order(agg3$dataset, agg3$scoring), ], row.names = FALSE)

  # Also count total XSS improvements per scoring mode (collapsed across datasets)
  message("\n--- Overall XSS improvement rate by scoring mode ---")
  overall <- aggregate(improved ~ scoring, data = xss, FUN = mean)
  names(overall)[2] <- "improvement_rate"
  print(overall, row.names = FALSE)

  invisible(list(xss_rates = agg, convergence = agg3, overall = overall))
}


# ---------------------------------------------------------------------------
# Run if called as script
# ---------------------------------------------------------------------------

if (!interactive() || identical(Sys.getenv("RUN_EXPERIMENT"), "true")) {
  message("Starting TBR <-> XSS cycling experiment...")
  t0 <- proc.time()
  results <- run_experiment(n_seeds = 20L)
  elapsed <- (proc.time() - t0)[3]
  message(sprintf("\nExperiment complete in %.1f seconds.", elapsed))

  saveRDS(results, "dev/benchmarks/expt_tbr_xss_results.rds")
  message("Results saved to dev/benchmarks/expt_tbr_xss_results.rds")

  analyse_results(results)
}
