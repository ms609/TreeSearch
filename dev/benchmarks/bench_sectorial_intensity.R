#!/usr/bin/env Rscript
# T-256: Sectorial search intensity experiment
#
# Compares four sectorial search configurations on gap datasets to determine
# if doubling/tripling XSS+RSS+CSS rounds closes the TNT score gap.
#
# DESIGNED FOR HAMILTON HPC. Do not run locally (60 runs × 30s ≈ 45 min).
#
# Usage:
#   Rscript bench_sectorial_intensity.R [timeout_s] [output_dir]
#
# Args:
#   timeout_s  — per-run budget in seconds (default: 30)
#   output_dir — where to write CSV results (default: .)

library(TreeSearch)
library(TreeTools)

args <- commandArgs(trailingOnly = TRUE)
timeout_s <- if (length(args) >= 1) as.integer(args[1]) else 30L
output_dir <- if (length(args) >= 2) args[2] else "."

cat("=== T-256: Sectorial Search Intensity Experiment ===\n")
cat(sprintf("Timeout: %ds\n", timeout_s))
cat(sprintf("TreeSearch version: %s\n", packageVersion("TreeSearch")))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

# ---- Datasets ----
# Five gap datasets from T-249/T-251 trajectory analysis
gap_names <- c("Conrad2008", "Geisler2001", "Wortley2006",
               "Zanol2014", "Zhu2013")

# Fitch mode: treat inapplicable as missing for TNT-fair comparison
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
# Each config overrides sectorial parameters from the "auto" preset.
# "auto" resolves to "default" for these 37–88 tip datasets (xss=3, rss=1, css=0).
configs <- list(
  baseline = list(
    label = "baseline",
    desc = "default preset as-is (xss=3, rss=1, css=0)"
  ),
  double = list(
    label = "2x_sect",
    desc = "doubled sectorial (xss=6, rss=2, css=1)",
    xssRounds = 6L, rssRounds = 2L, cssRounds = 1L
  ),
  triple = list(
    label = "3x_sect",
    desc = "tripled sectorial (xss=9, rss=3, css=2)",
    xssRounds = 9L, rssRounds = 3L, cssRounds = 2L
  ),
  nodrift_sect = list(
    label = "nodrift_3x",
    desc = "3x sectorial + no drift (xss=9, rss=3, css=2, drift=0)",
    xssRounds = 9L, rssRounds = 3L, cssRounds = 2L,
    driftCycles = 0L
  )
)

seeds <- c(1L, 2L, 3L)
total_runs <- length(configs) * length(datasets) * length(seeds)
cat(sprintf("Configs: %d, Datasets: %d, Seeds: %d → %d total runs\n\n",
            length(configs), length(datasets), length(seeds), total_runs))

# ---- TNT reference scores (from T-249 Fitch-mode benchmarks, 120s) ----
tnt_best <- c(
  Conrad2008 = 1725, Geisler2001 = 1293, Wortley2006 = 479,
  Zanol2014 = 1261, Zhu2013 = 624
)

# ---- Run experiments ----
results <- data.frame(
  dataset = character(), n_tips = integer(), n_chars = integer(),
  config = character(), seed = integer(), timeout_s = integer(),
  score = numeric(), n_trees = integer(), replicates = integer(),
  hits = integer(), wall_s = numeric(),
  tnt_best = numeric(), gap = numeric(),
  stringsAsFactors = FALSE
)

run_idx <- 0L
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  cat(sprintf("--- Config: %s (%s) ---\n", cfg$label, cfg$desc))

  for (ds_name in gap_names) {
    ds <- datasets[[ds_name]]
    ntip <- NTip(ds)
    nchar <- sum(attr(ds, "weight"))

    for (seed in seeds) {
      run_idx <- run_idx + 1L
      cat(sprintf("[%d/%d] %s / %s / seed=%d ... ",
                  run_idx, total_runs, cfg$label, ds_name, seed))

      set.seed(seed)

      # Build call args: strategy = "auto" + config-specific overrides
      call_args <- list(
        dataset = ds,
        concavity = Inf,
        maxReplicates = 96L,
        targetHits = 5L,
        maxSeconds = as.double(timeout_s),
        strategy = "auto",
        verbosity = 0L,
        nThreads = 1L
      )
      # Add sectorial parameter overrides (skip label/desc)
      override_names <- setdiff(names(cfg), c("label", "desc"))
      for (nm in override_names) {
        call_args[[nm]] <- cfg[[nm]]
      }

      t0 <- proc.time()
      result <- tryCatch(
        do.call(MaximizeParsimony, call_args),
        error = function(e) {
          warning("Error: ", ds_name, "/", cfg$label, ": ", conditionMessage(e))
          structure(list(), class = "multiPhylo",
                    score = NA_real_, pool_size = NA_integer_,
                    replicates = NA_integer_, hits_to_best = NA_integer_)
        }
      )
      wall_s <- as.double((proc.time() - t0)[3])

      sc <- attr(result, "score")
      tnt_ref <- tnt_best[ds_name]
      gap <- if (!is.na(sc)) sc - tnt_ref else NA_real_

      cat(sprintf("score=%s (gap=%s) in %.1fs (%d reps)\n",
                  if (is.na(sc)) "NA" else format(sc, nsmall = 0),
                  if (is.na(gap)) "NA" else sprintf("%+d", gap),
                  wall_s,
                  if (is.na(attr(result, "replicates"))) 0L
                  else attr(result, "replicates")))

      results <- rbind(results, data.frame(
        dataset = ds_name, n_tips = ntip, n_chars = nchar,
        config = cfg$label, seed = seed, timeout_s = timeout_s,
        score = sc, n_trees = length(result),
        replicates = if (is.na(attr(result, "replicates"))) NA_integer_
                     else attr(result, "replicates"),
        hits = if (is.na(attr(result, "hits_to_best"))) NA_integer_
               else attr(result, "hits_to_best"),
        wall_s = wall_s,
        tnt_best = tnt_ref, gap = gap,
        stringsAsFactors = FALSE
      ))
    }
  }
  cat("\n")
}

# ---- Write results ----
outfile <- file.path(output_dir,
                     sprintf("t256_sectorial_%ds_%s.csv",
                             timeout_s,
                             format(Sys.time(), "%Y%m%d_%H%M")))
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", outfile))

# ---- Summary ----
cat("\n=== Summary: Median score by config × dataset ===\n\n")
for (ds_name in gap_names) {
  cat(sprintf("  %s (TNT best: %d)\n", ds_name, tnt_best[ds_name]))
  for (cfg_name in names(configs)) {
    cfg <- configs[[cfg_name]]
    sub <- results[results$dataset == ds_name & results$config == cfg$label, ]
    med_score <- median(sub$score, na.rm = TRUE)
    med_gap <- median(sub$gap, na.rm = TRUE)
    best_score <- min(sub$score, na.rm = TRUE)
    med_reps <- median(sub$replicates, na.rm = TRUE)
    cat(sprintf("    %-15s median=%7.1f (gap %+5.1f)  best=%7.1f  reps=%.0f\n",
                cfg$label, med_score, med_gap, best_score, med_reps))
  }
  cat("\n")
}

# ---- Per-config aggregate ----
cat("=== Aggregate: Mean gap across all datasets ===\n\n")
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  sub <- results[results$config == cfg$label, ]
  mean_gap <- mean(sub$gap, na.rm = TRUE)
  med_gap <- median(sub$gap, na.rm = TRUE)
  cat(sprintf("  %-15s mean_gap=%+.1f  median_gap=%+.1f\n",
              cfg$label, mean_gap, med_gap))
}

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
