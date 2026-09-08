#!/usr/bin/env Rscript
# T-258: Intra-replicate fusing experiment
#
# Compares baseline vs intraFuse=TRUE on gap datasets to measure
# score quality and replicate throughput effects.
#
# DESIGNED FOR HAMILTON HPC. Do not run locally.
#
# Usage:
#   Rscript bench_intra_fuse.R [timeout_s] [output_dir]

library(TreeSearch)
library(TreeTools)

args <- commandArgs(trailingOnly = TRUE)
timeout_s <- if (length(args) >= 1) as.integer(args[1]) else 30L
output_dir <- if (length(args) >= 2) args[2] else "."

cat("=== T-258: Intra-Replicate Fusing Experiment ===\n")
cat(sprintf("Timeout: %ds\n", timeout_s))
cat(sprintf("TreeSearch version: %s\n", packageVersion("TreeSearch")))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

# ---- Datasets ----
gap_names <- c("Conrad2008", "Geisler2001", "Wortley2006",
               "Zanol2014", "Zhu2013")

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
  baseline = list(label = "baseline", desc = "default preset, no intra-fuse"),
  intra_fuse = list(label = "intra_fuse", desc = "default preset + intraFuse=TRUE",
                    intraFuse = TRUE)
)

seeds <- c(1L, 2L, 3L, 4L, 5L)  # 5 seeds for better signal
total_runs <- length(configs) * length(datasets) * length(seeds)
cat(sprintf("Configs: %d, Datasets: %d, Seeds: %d -> %d total runs\n\n",
            length(configs), length(datasets), length(seeds), total_runs))

# ---- TNT reference scores ----
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
                     sprintf("t258_intra_fuse_%ds_%s.csv",
                             timeout_s,
                             format(Sys.time(), "%Y%m%d_%H%M")))
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", outfile))

# ---- Summary ----
cat("\n=== Summary: Median score by config x dataset ===\n\n")
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

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
