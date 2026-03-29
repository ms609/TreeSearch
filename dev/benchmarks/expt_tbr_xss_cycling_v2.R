# Experiment v2: TBR <-> XSS cycling — EW vs IW, with timing + TAEB
#
# Question: After TBR converges, does XSS cycling pay for itself?
# Is the benefit stronger under IW than EW?
#
# Metric: Time-adjusted expected best (TAEB) — the expected minimum score
# from k = floor(budget / per_rep_time) independent replicates.
# This is the right metric because conditions A and B have different
# per-replicate costs.
#
# Protocol per (dataset x scoring x seed):
#   Shared: Random Wagner -> NNI -> TBR convergence
#   Condition A (TBR-only): stop here
#   Condition B (TBR+XSS): [XSS -> TBR] cycling until joint convergence
#
# Same seed -> same starting tree -> paired comparison.
#
# Usage:
#   Rscript dev/benchmarks/expt_tbr_xss_cycling_v2.R

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
  switch(mode,
    EW   = Inf,
    IW10 = 10,
    IW3  = 3,
    stop("Unknown scoring mode: ", mode)
  )
}

# TAEB: expected minimum from k independent draws
# scores: numeric vector of per-replicate scores
# k: number of replicates fitting in budget
# B: number of bootstrap resamples
taeb <- function(scores, k, B = 5000L) {
  if (k < 1) return(NA_real_)
  if (k >= length(scores)) return(min(scores))
  mins <- replicate(B, min(sample(scores, k, replace = TRUE)))
  mean(mins)
}

# Run one seed: shared baseline + both conditions
run_one <- function(ds_data, min_steps, scoring_mode, seed,
                    max_cycles = 10L) {
  set.seed(seed)
  k <- concavity_value(scoring_mode)

  # --- Shared baseline: Wagner -> NNI -> TBR ---
  t0 <- proc.time()[3]

  wagner <- TreeSearch:::ts_random_wagner_tree(
    ds_data$contrast, ds_data$tip_data, ds_data$weight, ds_data$levels,
    min_steps = min_steps, concavity = k
  )

  nni <- TreeSearch:::ts_nni_search(
    wagner$edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
    ds_data$levels, maxHits = 20L, min_steps = min_steps, concavity = k
  )

  tbr <- TreeSearch:::ts_tbr_search(
    nni$edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
    ds_data$levels, maxHits = 1L, min_steps = min_steps, concavity = k
  )

  tbr_time <- proc.time()[3] - t0

  # Condition A: TBR-only result
  condA <- data.frame(
    condition    = "TBR_only",
    score        = tbr$score,
    time_s       = tbr_time,
    n_xss_cycles = 0L,
    xss_improved = 0L,
    stringsAsFactors = FALSE
  )

  # --- Condition B: XSS <-> TBR cycling from TBR-converged tree ---
  t1 <- proc.time()[3]
  current_edge  <- tbr$edge
  current_score <- tbr$score
  n_xss_cycles  <- 0L
  xss_improved  <- 0L

  for (cyc in seq_len(max_cycles)) {
    xss <- TreeSearch:::ts_xss_search(
      current_edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
      ds_data$levels,
      nPartitions = 4L, xssRounds = 3L, acceptEqual = FALSE,
      maxHits = 1L, min_steps = min_steps, concavity = k
    )

    xss_got_better <- xss$score < current_score - 1e-9
    n_xss_cycles <- n_xss_cycles + 1L
    if (xss_got_better) xss_improved <- xss_improved + 1L

    current_edge  <- xss$edge
    current_score <- xss$score

    tbr2 <- TreeSearch:::ts_tbr_search(
      current_edge, ds_data$contrast, ds_data$tip_data, ds_data$weight,
      ds_data$levels, maxHits = 1L, min_steps = min_steps, concavity = k
    )

    tbr_got_better <- tbr2$score < current_score - 1e-9

    current_edge  <- tbr2$edge
    current_score <- tbr2$score

    if (!xss_got_better && !tbr_got_better) break
  }

  cycling_time <- proc.time()[3] - t1

  condB <- data.frame(
    condition    = "TBR_XSS",
    score        = current_score,
    time_s       = tbr_time + cycling_time,
    n_xss_cycles = n_xss_cycles,
    xss_improved = xss_improved,
    stringsAsFactors = FALSE
  )

  rbind(condA, condB)
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

run_experiment <- function(n_seeds = 20L) {

  # Datasets: 4 standard + 1 large
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
  message("  TBR <-> XSS Cycling: Definitive Results")
  message("============================================================")

  # --- 1. XSS improvement rate ---
  message("\n--- 1. XSS improvement rate (condition B, by dataset x scoring) ---")
  condB <- df[df$condition == "TBR_XSS", ]
  condB$improved <- condB$xss_improved > 0
  agg_rate <- aggregate(improved ~ dataset + scoring + n_tip,
                        data = condB, FUN = mean)
  agg_rate$rate_pct <- round(agg_rate$improved * 100, 1)
  print(agg_rate[order(agg_rate$n_tip, agg_rate$scoring),
                 c("dataset", "n_tip", "scoring", "rate_pct")],
        row.names = FALSE)

  # --- 2. Score delta (B - A, paired by seed) ---
  message("\n--- 2. Paired score delta (TBR_XSS - TBR_only) ---")
  condA <- df[df$condition == "TBR_only", ]
  condB <- df[df$condition == "TBR_XSS", ]
  paired <- merge(
    condA[, c("dataset", "scoring", "seed", "score", "time_s")],
    condB[, c("dataset", "scoring", "seed", "score", "time_s")],
    by = c("dataset", "scoring", "seed"),
    suffixes = c("_A", "_B")
  )
  paired$delta_score <- paired$score_B - paired$score_A
  paired$delta_time  <- paired$time_s_B - paired$time_s_A
  paired$overhead_pct <- round(100 * paired$delta_time / paired$time_s_A, 1)

  agg_delta <- aggregate(
    cbind(delta_score, overhead_pct) ~ dataset + scoring,
    data = paired, FUN = function(x) round(mean(x), 2)
  )
  agg_delta <- agg_delta[order(agg_delta$dataset, agg_delta$scoring), ]
  names(agg_delta)[3:4] <- c("mean_delta_score", "mean_overhead_pct")
  print(agg_delta, row.names = FALSE)

  # --- 3. Time-adjusted expected best (TAEB) ---
  message("\n--- 3. TAEB at multiple budgets ---")
  budgets <- c(10, 30, 60, 120)
  taeb_rows <- list()
  ri <- 0L

  combos <- unique(df[, c("dataset", "scoring")])
  for (i in seq_len(nrow(combos))) {
    ds_name <- combos$dataset[i]
    sc_mode <- combos$scoring[i]

    for (cond in c("TBR_only", "TBR_XSS")) {
      sub <- df[df$dataset == ds_name & df$scoring == sc_mode &
                  df$condition == cond, ]
      scores <- sub$score
      med_time <- median(sub$time_s)

      for (b in budgets) {
        k <- floor(b / med_time)
        ri <- ri + 1L
        taeb_rows[[ri]] <- data.frame(
          dataset    = ds_name,
          scoring    = sc_mode,
          condition  = cond,
          budget_s   = b,
          med_time_s = round(med_time, 2),
          k_reps     = k,
          taeb       = round(taeb(scores, k), 2),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  taeb_df <- do.call(rbind, taeb_rows)

  # Pivot: show A vs B side by side per (dataset, scoring, budget)
  taeb_A <- taeb_df[taeb_df$condition == "TBR_only",
                     c("dataset", "scoring", "budget_s",
                       "med_time_s", "k_reps", "taeb")]
  taeb_B <- taeb_df[taeb_df$condition == "TBR_XSS",
                     c("dataset", "scoring", "budget_s",
                       "med_time_s", "k_reps", "taeb")]
  names(taeb_A)[4:6] <- paste0(names(taeb_A)[4:6], "_A")
  names(taeb_B)[4:6] <- paste0(names(taeb_B)[4:6], "_B")

  taeb_wide <- merge(taeb_A, taeb_B,
                     by = c("dataset", "scoring", "budget_s"))
  taeb_wide$delta_taeb <- taeb_wide$taeb_B - taeb_wide$taeb_A

  # Negative delta = B is better (lower score wins)
  taeb_wide <- taeb_wide[order(taeb_wide$dataset,
                               taeb_wide$scoring,
                               taeb_wide$budget_s), ]

  message("\nFull TAEB table (delta < 0 means XSS cycling helps):")
  print(taeb_wide[, c("dataset", "scoring", "budget_s",
                       "k_reps_A", "taeb_A", "k_reps_B", "taeb_B",
                       "delta_taeb")],
        row.names = FALSE)

  # --- 4. Summary verdict per scoring mode ---
  message("\n--- 4. Summary: mean TAEB delta by scoring mode + budget ---")
  agg_verdict <- aggregate(delta_taeb ~ scoring + budget_s,
                           data = taeb_wide, FUN = mean)
  agg_verdict <- agg_verdict[order(agg_verdict$scoring,
                                   agg_verdict$budget_s), ]
  names(agg_verdict)[3] <- "mean_delta_taeb"
  agg_verdict$mean_delta_taeb <- round(agg_verdict$mean_delta_taeb, 2)
  print(agg_verdict, row.names = FALSE)

  # --- 5. Per-dataset verdict (averaged over budgets) ---
  message("\n--- 5. Mean TAEB delta by dataset x scoring (across budgets) ---")
  agg_ds <- aggregate(delta_taeb ~ dataset + scoring,
                      data = taeb_wide, FUN = mean)
  agg_ds <- agg_ds[order(agg_ds$dataset, agg_ds$scoring), ]
  agg_ds$delta_taeb <- round(agg_ds$delta_taeb, 2)
  print(agg_ds, row.names = FALSE)

  invisible(list(
    raw = df,
    paired = paired,
    taeb = taeb_wide,
    verdict_mode = agg_verdict,
    verdict_ds = agg_ds
  ))
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if (!interactive() || identical(Sys.getenv("RUN_EXPERIMENT"), "true")) {
  message("Starting TBR <-> XSS cycling experiment v2...")
  t0 <- proc.time()
  results <- run_experiment(n_seeds = 20L)
  elapsed <- (proc.time() - t0)[3]
  message(sprintf("\nExperiment complete in %.1f seconds.", elapsed))

  saveRDS(results, "dev/benchmarks/expt_tbr_xss_v2_results.rds")
  message("Results saved to dev/benchmarks/expt_tbr_xss_v2_results.rds")

  out <- analyse_results(results)
}
