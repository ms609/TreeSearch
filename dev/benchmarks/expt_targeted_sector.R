# Experiment: Targeted post-clip sectorial search
#
# After each accepted TBR move, run a sector-masked TBR search on the
# just-moved clip subtree. Does refining the subtree's internal
# arrangement in its new context improve scores? Is the effect stronger
# under IW?
#
# Condition A: Wagner -> NNI -> TBR (normal)
# Condition B: Wagner -> NNI -> TBR (with targetedSector = TRUE)
# Same seed -> same starting tree -> paired comparison.

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
  as.integer(MinimumLength(ds, compress = TRUE))
}

concavity_value <- function(mode) {
  switch(mode, EW = Inf, IW10 = 10, IW3 = 3, stop("Unknown: ", mode))
}

taeb <- function(scores, k, B = 5000L) {
  if (k < 1) return(NA_real_)
  if (k >= length(scores)) return(min(scores))
  mins <- replicate(B, min(sample(scores, k, replace = TRUE)))
  mean(mins)
}

# Run one seed: shared Wagner+NNI start, then A (normal TBR) and B (targeted TBR)
run_one <- function(ds_data, min_steps, scoring_mode, seed) {
  set.seed(seed)
  k <- concavity_value(scoring_mode)

  wagner <- TreeSearch:::ts_random_wagner_tree(
    ds_data$contrast, ds_data$tip_data, ds_data$weight, ds_data$levels,
    min_steps = min_steps, concavity = k
  )

  nni <- TreeSearch:::ts_nni_search(
    wagner$edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
    ds_data$levels, maxHits = 20L, min_steps = min_steps, concavity = k
  )

  nni_edge <- nni$edge

  # --- Condition A: normal TBR ---
  t0 <- proc.time()[3]
  tbrA <- TreeSearch:::ts_tbr_search(
    nni_edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
    ds_data$levels, maxHits = 1L, min_steps = min_steps, concavity = k,
    targetedSector = FALSE
  )
  time_A <- proc.time()[3] - t0

  # --- Condition B: targeted TBR ---
  # Re-seed to get the same NNI starting tree
  # Actually, we use the same nni_edge for both conditions
  t1 <- proc.time()[3]
  tbrB <- TreeSearch:::ts_tbr_search(
    nni_edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
    ds_data$levels, maxHits = 1L, min_steps = min_steps, concavity = k,
    targetedSector = TRUE, targetedMinSize = 6L
  )
  time_B <- proc.time()[3] - t1

  rbind(
    data.frame(
      condition = "normal",
      score = tbrA$score,
      time_s = time_A,
      n_accepted = tbrA$n_accepted,
      n_targeted_calls = tbrA$n_targeted_calls,
      n_targeted_improved = tbrA$n_targeted_improved,
      stringsAsFactors = FALSE
    ),
    data.frame(
      condition = "targeted",
      score = tbrB$score,
      time_s = time_B,
      n_accepted = tbrB$n_accepted,
      n_targeted_calls = tbrB$n_targeted_calls,
      n_targeted_improved = tbrB$n_targeted_improved,
      stringsAsFactors = FALSE
    )
  )
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

run_experiment <- function(n_seeds = 20L) {

  ds_specs <- list(
    list(name = "Agnarsson2004", src = "builtin"),
    list(name = "Zhu2013",       src = "builtin"),
    list(name = "Giles2015",     src = "builtin"),
    list(name = "Dikow2009",     src = "builtin"),
    list(name = "mbank_X30754",  src = "dev/benchmarks/mbank_X30754.nex")
  )

  scoring_modes <- c("EW", "IW10", "IW3")
  seeds <- sample.int(10000, n_seeds)

  all_results <- vector("list",
    length(ds_specs) * length(scoring_modes) * n_seeds)
  idx <- 0L

  for (spec in ds_specs) {
    ds_name <- spec$name
    message("=== ", ds_name, " ===")

    if (spec$src == "builtin") {
      dataset <- TreeSearch::inapplicable.phyData[[ds_name]]
    } else {
      dataset <- ReadAsPhyDat(spec$src)
    }
    ds_data   <- make_ts_data(dataset)
    min_steps <- min_steps_for(dataset)
    n_tip     <- length(dataset)

    for (mode in scoring_modes) {
      message("  ", mode, ": ", appendLF = FALSE)
      for (s in seq_along(seeds)) {
        if (s %% 5 == 0) message(s, " ", appendLF = FALSE)
        idx <- idx + 1L
        res <- run_one(ds_data, min_steps, mode, seeds[s])
        res$dataset <- ds_name
        res$scoring <- mode
        res$seed    <- seeds[s]
        res$n_tip   <- n_tip
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

  message("\n============================================================")
  message("  Targeted Post-Clip Sector Search: Results")
  message("============================================================")

  condB <- df[df$condition == "targeted", ]

  # --- 1. Targeting diagnostics ---
  message("\n--- 1. Targeting diagnostics (condition B only) ---")
  agg_diag <- aggregate(
    cbind(n_accepted, n_targeted_calls, n_targeted_improved) ~
      dataset + scoring + n_tip,
    data = condB, FUN = mean
  )
  agg_diag$hit_rate <- ifelse(agg_diag$n_targeted_calls > 0,
    round(agg_diag$n_targeted_improved / agg_diag$n_targeted_calls * 100, 1),
    0)
  agg_diag$pct_eligible <- round(
    agg_diag$n_targeted_calls / agg_diag$n_accepted * 100, 1)
  agg_diag <- agg_diag[order(agg_diag$n_tip, agg_diag$scoring), ]
  print(agg_diag[, c("dataset", "n_tip", "scoring", "n_accepted",
                      "n_targeted_calls", "n_targeted_improved",
                      "hit_rate", "pct_eligible")],
        row.names = FALSE)

  # --- 2. Paired score delta ---
  message("\n--- 2. Paired score delta (targeted - normal) ---")
  condA <- df[df$condition == "normal", ]
  paired <- merge(
    condA[, c("dataset", "scoring", "seed", "score", "time_s")],
    condB[, c("dataset", "scoring", "seed", "score", "time_s")],
    by = c("dataset", "scoring", "seed"),
    suffixes = c("_A", "_B")
  )
  paired$delta_score <- paired$score_B - paired$score_A
  paired$overhead_pct <- round(100 * (paired$time_s_B - paired$time_s_A) /
                                 paired$time_s_A, 1)

  agg_delta <- aggregate(
    cbind(delta_score, overhead_pct) ~ dataset + scoring,
    data = paired, FUN = function(x) round(mean(x), 3)
  )
  agg_delta <- agg_delta[order(agg_delta$dataset, agg_delta$scoring), ]
  names(agg_delta)[3:4] <- c("mean_delta", "mean_overhead_pct")
  print(agg_delta, row.names = FALSE)

  # Count how many seeds improved
  message("\n--- Seeds improved (out of 20) ---")
  agg_improved <- aggregate(delta_score ~ dataset + scoring, data = paired,
    FUN = function(x) sum(x < -1e-9))
  names(agg_improved)[3] <- "n_improved"
  agg_improved <- agg_improved[order(agg_improved$dataset,
                                     agg_improved$scoring), ]
  print(agg_improved, row.names = FALSE)

  # --- 3. TAEB ---
  message("\n--- 3. TAEB at multiple budgets ---")
  budgets <- c(10, 30, 60, 120)
  taeb_rows <- list()
  ri <- 0L

  combos <- unique(df[, c("dataset", "scoring")])
  for (i in seq_len(nrow(combos))) {
    ds_name <- combos$dataset[i]
    sc_mode <- combos$scoring[i]

    for (cond in c("normal", "targeted")) {
      sub <- df[df$dataset == ds_name & df$scoring == sc_mode &
                  df$condition == cond, ]
      med_time <- median(sub$time_s)

      for (b in budgets) {
        k <- floor(b / med_time)
        ri <- ri + 1L
        taeb_rows[[ri]] <- data.frame(
          dataset = ds_name, scoring = sc_mode, condition = cond,
          budget_s = b, med_time_s = round(med_time, 2),
          k_reps = k, taeb = round(taeb(sub$score, k), 2),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  taeb_df <- do.call(rbind, taeb_rows)

  taeb_A <- taeb_df[taeb_df$condition == "normal",
                     c("dataset", "scoring", "budget_s", "k_reps", "taeb")]
  taeb_B <- taeb_df[taeb_df$condition == "targeted",
                     c("dataset", "scoring", "budget_s", "k_reps", "taeb")]
  names(taeb_A)[4:5] <- paste0(names(taeb_A)[4:5], "_A")
  names(taeb_B)[4:5] <- paste0(names(taeb_B)[4:5], "_B")

  taeb_wide <- merge(taeb_A, taeb_B, by = c("dataset", "scoring", "budget_s"))
  taeb_wide$delta_taeb <- taeb_wide$taeb_B - taeb_wide$taeb_A
  taeb_wide <- taeb_wide[order(taeb_wide$dataset, taeb_wide$scoring,
                               taeb_wide$budget_s), ]

  message("\nFull TAEB table (delta < 0 means targeting helps):")
  print(taeb_wide[, c("dataset", "scoring", "budget_s",
                       "k_reps_A", "taeb_A", "k_reps_B", "taeb_B",
                       "delta_taeb")],
        row.names = FALSE)

  # --- 4. Summary by scoring mode ---
  message("\n--- 4. Mean TAEB delta by scoring mode + budget ---")
  agg_verdict <- aggregate(delta_taeb ~ scoring + budget_s,
                           data = taeb_wide, FUN = mean)
  agg_verdict$delta_taeb <- round(agg_verdict$delta_taeb, 2)
  agg_verdict <- agg_verdict[order(agg_verdict$scoring,
                                   agg_verdict$budget_s), ]
  print(agg_verdict, row.names = FALSE)

  # --- 5. IW vs EW targeting hit rate comparison ---
  message("\n--- 5. Targeting hit rate: IW vs EW ---")
  hr <- aggregate(
    cbind(n_targeted_calls, n_targeted_improved) ~ scoring,
    data = condB, FUN = sum
  )
  hr$hit_rate_pct <- round(hr$n_targeted_improved / hr$n_targeted_calls * 100, 1)
  print(hr, row.names = FALSE)

  invisible(list(raw = df, paired = paired, taeb = taeb_wide,
                 diagnostics = agg_diag))
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if (!interactive() || identical(Sys.getenv("RUN_EXPERIMENT"), "true")) {
  message("Starting targeted post-clip sector search experiment...")
  t0 <- proc.time()
  results <- run_experiment(n_seeds = 20L)
  elapsed <- (proc.time() - t0)[3]
  message(sprintf("\nExperiment complete in %.1f seconds.", elapsed))

  saveRDS(results, "dev/benchmarks/expt_targeted_sector_results.rds")
  message("Results saved to dev/benchmarks/expt_targeted_sector_results.rds")

  out <- analyse_results(results)
}
