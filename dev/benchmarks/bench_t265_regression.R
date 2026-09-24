#!/usr/bin/env Rscript
# T-265: Per-replicate search quality regression diagnosis
#
# DESIGNED FOR HAMILTON HPC. Do not run locally.
#
# Tests whether the quality regression is in preset params vs engine code
# by comparing 3 configurations on the datasets with largest TNT gaps:
#
#   r2_equiv     — Minimal pipeline matching R2 structure: 12 ratchet (4%,
#                  auto moves), 2 drift, no sectorial, 1 Wagner, no tabu,
#                  no NNI warmup. Tests what R2 actually ran.
#   r2_modern    — R2 structure + modern ratchet tuning: 12 ratchet (25%,
#                  5 moves), 0 drift, 1 Wagner, no sectorial, no tabu,
#                  NNI warmup ON. Tests whether modern ratchet params help
#                  with a minimal pipeline.
#   auto_preset  — Current auto-selected preset (default or thorough).
#                  Tests whether added complexity helps or hurts.
#
# If r2_equiv or r2_modern produce better scores -> preset complexity is
# the problem. If all configs show the same regression -> engine code issue.
#
# Usage:
#   Rscript bench_t265_regression.R [timeout_s] [output_dir]
#
# Default: 120s budget, output to current directory.

library(TreeSearch)
library(TreeTools)

args <- commandArgs(trailingOnly = TRUE)
timeout_s <- if (length(args) >= 1) as.integer(args[1]) else 120L
output_dir <- if (length(args) >= 2) args[2] else "."

cat("=== T-265: Per-Replicate Quality Regression Diagnosis ===\n")
cat(sprintf("Timeout: %ds\n", timeout_s))
cat(sprintf("TreeSearch version: %s\n", packageVersion("TreeSearch")))
cat(sprintf("Output dir: %s\n", output_dir))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

# ---- Datasets ----
# 8 datasets with largest persistent TNT gaps, plus Wilson2003 from T-265
gap_names <- c(
  "Wortley2006", "Eklund2004", "Wilson2003", "Conrad2008",
  "Geisler2001", "Zanol2014", "Zhu2013", "Giles2015", "Dikow2009"
)

# Convert inapplicable to missing for EW Fitch scoring (match TNT)
fitch_mode <- function(dataset) {
  contrast <- attr(dataset, "contrast")
  levels <- attr(dataset, "levels")
  inapp_col <- match("-", levels)
  if (is.na(inapp_col)) return(dataset)
  for (i in seq_len(nrow(contrast))) {
    if (contrast[i, inapp_col] == 1 && sum(contrast[i, ]) == 1) {
      contrast[i, ] <- 1
    }
  }
  attr(dataset, "contrast") <- contrast
  dataset
}

datasets <- lapply(
  setNames(gap_names, gap_names),
  function(nm) fitch_mode(inapplicable.phyData[[nm]])
)

# ---- Configurations ----
configs <- list(
  r2_equiv = list(
    label = "r2_equiv",
    desc = "R2 pipeline: 12 ratchet (4%), 2 drift, no sectorial, no tabu",
    control = SearchControl(
      ratchetCycles = 12L,
      ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L,
      ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 2L,
      driftAfdLimit = 5L,
      driftRfdLimit = 0.15,
      xssRounds = 0L, rssRounds = 0L, cssRounds = 0L,
      wagnerStarts = 1L,
      tabuSize = 0L,
      nniFirst = FALSE, sprFirst = FALSE,
      perturbStopFactor = 0L,
      adaptiveLevel = FALSE,
      maxOuterResets = 0L,
      outerCycles = 1L,
      fuseInterval = 5L,
      fuseAcceptEqual = FALSE,
      poolMaxSize = 100L,
      consensusStableReps = 0L,
      nniPerturbCycles = 0L,
      annealCycles = 0L,
      adaptiveStart = FALSE,
      enumTimeFraction = 0.1
    )
  ),
  r2_modern = list(
    label = "r2_modern",
    desc = "R2 structure + modern ratchet (25%, 5 moves), NNI warmup, no drift",
    control = SearchControl(
      ratchetCycles = 12L,
      ratchetPerturbProb = 0.25,
      ratchetPerturbMode = 0L,
      ratchetPerturbMaxMoves = 5L,
      ratchetAdaptive = FALSE,
      driftCycles = 0L,
      xssRounds = 0L, rssRounds = 0L, cssRounds = 0L,
      wagnerStarts = 1L,
      tabuSize = 0L,
      nniFirst = TRUE, sprFirst = FALSE,
      perturbStopFactor = 0L,
      adaptiveLevel = FALSE,
      maxOuterResets = 0L,
      outerCycles = 1L,
      fuseInterval = 5L,
      fuseAcceptEqual = FALSE,
      poolMaxSize = 100L,
      consensusStableReps = 0L,
      nniPerturbCycles = 0L,
      annealCycles = 0L,
      adaptiveStart = FALSE,
      enumTimeFraction = 0.1
    )
  ),
  auto_preset = list(
    label = "auto_preset",
    desc = "Current auto-selected preset (default or thorough)"
    # No control override — uses strategy = "auto"
  )
)

seeds <- 1:5
total_runs <- length(configs) * length(datasets) * length(seeds)
cat(sprintf("Configs: %d, Datasets: %d, Seeds: %d -> %d total runs\n",
            length(configs), length(datasets), length(seeds), total_runs))

# TNT reference scores (from bench_intra_fuse.R and T-265 notes)
tnt_best <- c(
  Wortley2006 = 479, Eklund2004 = 438, Wilson2003 = 860,
  Conrad2008 = 1725, Geisler2001 = 1293,
  Zanol2014 = 1261, Zhu2013 = 624,
  Giles2015 = 670, Dikow2009 = 1603
)

# ---- Run experiments ----
results <- data.frame(
  dataset = character(), n_tips = integer(), n_patterns = integer(),
  auto_strategy = character(),
  config = character(), seed = integer(), timeout_s = integer(),
  score = numeric(), n_trees = integer(), replicates = integer(),
  hits = integer(), wall_s = numeric(),
  tnt_best = numeric(), gap = numeric(),
  stringsAsFactors = FALSE
)

run_idx <- 0L
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  cat(sprintf("\n--- Config: %s (%s) ---\n", cfg$label, cfg$desc))

  for (ds_name in gap_names) {
    ds <- datasets[[ds_name]]
    ntip <- length(ds)
    npat <- sum(attr(ds, "weight"))
    auto_strat <- if (ntip <= 30) "sprint"
      else if (npat < 100) "default"
      else if (ntip >= 120) "large"
      else if (ntip >= 65) "thorough"
      else "default"

    for (s in seeds) {
      run_idx <- run_idx + 1L
      cat(sprintf("  [%d/%d] %s / %s / seed=%d ... ",
                  run_idx, total_runs, ds_name, cfg$label, s))

      set.seed(s)
      t0 <- proc.time()

      tryCatch({
        if (cfg_name == "auto_preset") {
          res <- MaximizeParsimony(
            ds,
            maxSeconds = timeout_s,
            strategy = "auto",
            verbosity = 0L,
            nThreads = 1L
          )
        } else {
          res <- MaximizeParsimony(
            ds,
            maxSeconds = timeout_s,
            strategy = "none",
            control = cfg$control,
            verbosity = 0L,
            nThreads = 1L
          )
        }

        elapsed <- (proc.time() - t0)[3]
        best_score <- attr(res, "score")
        n_trees <- length(res)
        reps <- attr(res, "replicates")
        hits <- attr(res, "hits")
        tnt_ref <- tnt_best[ds_name]
        gap <- if (!is.na(tnt_ref)) best_score - tnt_ref else NA_real_

        cat(sprintf("score=%g, gap=%s, reps=%d, %.1fs\n",
                    best_score,
                    if (is.na(gap)) "?" else sprintf("%+d", gap),
                    reps, elapsed))

        results <- rbind(results, data.frame(
          dataset = ds_name, n_tips = ntip, n_patterns = npat,
          auto_strategy = auto_strat,
          config = cfg$label, seed = s, timeout_s = timeout_s,
          score = best_score, n_trees = n_trees, replicates = reps,
          hits = hits, wall_s = elapsed,
          tnt_best = tnt_ref, gap = gap,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        elapsed <- (proc.time() - t0)[3]
        cat(sprintf("ERROR: %s (%.1fs)\n", conditionMessage(e), elapsed))
        results <<- rbind(results, data.frame(
          dataset = ds_name, n_tips = ntip, n_patterns = npat,
          auto_strategy = auto_strat,
          config = cfg$label, seed = s, timeout_s = timeout_s,
          score = NA_real_, n_trees = NA_integer_, replicates = NA_integer_,
          hits = NA_integer_, wall_s = elapsed,
          tnt_best = tnt_best[ds_name], gap = NA_real_,
          stringsAsFactors = FALSE
        ))
      })
    }
  }
}

# ---- Save results ----
out_file <- file.path(output_dir,
                      sprintf("t265_results_%ds.csv", timeout_s))
write.csv(results, out_file, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", out_file))

# ---- Summary ----
cat("\n=== Summary by config × dataset (median score, median gap) ===\n\n")
for (ds_name in gap_names) {
  sub <- results[results$dataset == ds_name, ]
  if (nrow(sub) == 0) next
  cat(sprintf("  %s (%dt, %dp, auto=%s, TNT=%s):\n",
              ds_name, sub$n_tips[1], sub$n_patterns[1],
              sub$auto_strategy[1],
              if (is.na(tnt_best[ds_name])) "?" else tnt_best[ds_name]))
  for (cfg_name in names(configs)) {
    cfg_sub <- sub[sub$config == configs[[cfg_name]]$label, ]
    if (nrow(cfg_sub) == 0) next
    med_score <- median(cfg_sub$score, na.rm = TRUE)
    med_gap <- median(cfg_sub$gap, na.rm = TRUE)
    min_score <- min(cfg_sub$score, na.rm = TRUE)
    max_score <- max(cfg_sub$score, na.rm = TRUE)
    med_reps <- median(cfg_sub$replicates, na.rm = TRUE)
    unique_scores <- length(unique(na.omit(cfg_sub$score)))
    cat(sprintf("    %-14s median=%7.0f (range %g-%g), gap=%+.0f, reps=%.0f, unique_scores=%d\n",
                configs[[cfg_name]]$label, med_score, min_score, max_score,
                med_gap, med_reps, unique_scores))
  }
}

# ---- Per-replicate convergence check ----
cat("\n=== Score diversity across seeds (do all seeds find the same score?) ===\n\n")
for (ds_name in gap_names) {
  sub <- results[results$dataset == ds_name, ]
  if (nrow(sub) == 0) next
  cat(sprintf("  %s:\n", ds_name))
  for (cfg_name in names(configs)) {
    cfg_sub <- sub[sub$config == configs[[cfg_name]]$label, ]
    if (nrow(cfg_sub) == 0) next
    scores <- na.omit(cfg_sub$score)
    if (length(scores) == 0) next
    n_unique <- length(unique(scores))
    cat(sprintf("    %-14s scores: %s  (%d unique)\n",
                configs[[cfg_name]]$label,
                paste(scores, collapse = ", "),
                n_unique))
  }
}

cat(sprintf("\nCompleted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
