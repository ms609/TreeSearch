#!/usr/bin/env Rscript
# T-289: Prune-reinsert perturbation benchmark
#
# DESIGNED FOR HAMILTON HPC. Do not run locally (hours of wall time).
#
# Evaluates taxon prune-reinsert (T-266) as perturbation strategy:
#   - Does it improve scores vs baseline (no prune-reinsert)?
#   - Optimal cycle count (1, 3, 5)?
#   - Optimal drop fraction (0.05, 0.10, 0.20, 0.30)?
#   - RANDOM vs INSTABILITY tip selection?
#   - Complement or substitute for ratchet?
#   - Scaling with dataset size (37t → 180t)?
#
# Design: Two-stage grid.
#   Stage 1 ("sweep"): coarse grid on 5 datasets, 5 seeds, 30s budget.
#     Identifies best cycle count and drop fraction.
#   Stage 2 ("confirm"): best configs + baseline on 5 datasets,
#     5 seeds, 30s + 60s budgets. Also tests INSTABILITY selection
#     and ratchet-replacement.
#
# Usage:
#   Rscript bench_prune_reinsert.R [stage] [timeout_s] [output_dir]
#   stage: 1 (sweep) or 2 (confirm). Default: 1
#   timeout_s: search budget in seconds. Default: 30
#   output_dir: where to write CSV results. Default: "."
#
# Output: t289_stage{1,2}_{timeout}s.csv

library(TreeSearch)
library(TreeTools)

args <- commandArgs(trailingOnly = TRUE)
stage     <- if (length(args) >= 1) as.integer(args[1]) else 1L
timeout_s <- if (length(args) >= 2) as.integer(args[2]) else 30L
output_dir <- if (length(args) >= 3) args[3] else "."

cat("=== T-289: Prune-Reinsert Benchmark ===\n")
cat(sprintf("Stage: %d, Timeout: %ds\n", stage, timeout_s))
cat(sprintf("TreeSearch version: %s\n", packageVersion("TreeSearch")))
cat(sprintf("Output dir: %s\n", output_dir))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

# ---- Datasets ----
# 5 datasets spanning tip-count range, chosen for enough landscape
# difficulty that perturbation quality matters.
bench_names <- c(
  "Wortley2006",    #  37 tips — small, gap dataset
  "Agnarsson2004",  #  62 tips — medium
  "Zhu2013",        #  75 tips — hard, high missing
  "Dikow2009"       #  88 tips — largest standard
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
  setNames(bench_names, bench_names),
  function(nm) fitch_mode(inapplicable.phyData[[nm]])
)

# Also load 180-tip dataset if available
mbank_path <- file.path(dirname(dirname(getwd())),
                        "TreeSearch-a", "dev", "benchmarks",
                        "mbank_X30754.nex")
if (!file.exists(mbank_path)) {
  mbank_path <- Sys.glob("/nobackup/*/TreeSearch-a/dev/benchmarks/mbank_X30754.nex")
  if (length(mbank_path) > 0) mbank_path <- mbank_path[1]
}
if (length(mbank_path) == 1 && file.exists(mbank_path)) {
  cat("Loading 180-tip dataset from:", mbank_path, "\n")
  mbank <- fitch_mode(ReadAsPhyDat(mbank_path))
  datasets[["mbank_X30754"]] <- mbank
  bench_names <- c(bench_names, "mbank_X30754")
} else {
  cat("180-tip dataset not found; skipping mbank_X30754\n")
}

# TNT reference scores (EW Fitch mode)
tnt_best <- c(
  Wortley2006 = 479, Agnarsson2004 = 718,
  Zhu2013 = 624, Dikow2009 = 1603,
  mbank_X30754 = 1164
)

seeds <- 1:5

# ---- Baseline control (current auto preset, no prune-reinsert) ----
make_baseline <- function() {
  # No prune-reinsert; everything else at default
  SearchControl(
    pruneReinsertCycles = 0L,
    consensusStableReps = 0L,
    nniPerturbCycles = 0L,
    driftCycles = 0L
  )
}

# ---- Build config grid ----
build_configs <- function(stage) {
  cfgs <- list()

  # Baseline: no prune-reinsert
  cfgs[["baseline"]] <- list(
    label = "baseline",
    desc = "No prune-reinsert (auto preset)",
    control = NULL  # use strategy = "auto"
  )

  if (stage == 1L) {
    # Stage 1: sweep cycles × drop_fraction, RANDOM selection only
    for (cyc in c(1L, 3L, 5L)) {
      for (drop in c(0.05, 0.10, 0.20, 0.30)) {
        tag <- sprintf("pr_c%d_d%02d", cyc, round(drop * 100))
        cfgs[[tag]] <- list(
          label = tag,
          desc = sprintf("PR: %d cycles, %.0f%% drop, random",
                         cyc, drop * 100),
          pr_cycles = cyc,
          pr_drop = drop,
          pr_selection = 0L  # RANDOM
        )
      }
    }
  } else {
    # Stage 2: best configs from stage 1 + INSTABILITY + ratchet-replacement
    # (Placeholder — human fills in best configs after stage 1 analysis)
    # Default stage 2 tests a reasonable spread:
    for (cyc in c(1L, 3L)) {
      for (drop in c(0.10, 0.20)) {
        for (sel in c(0L, 1L)) {
          sel_tag <- if (sel == 0L) "rand" else "inst"
          tag <- sprintf("pr_c%d_d%02d_%s", cyc, round(drop * 100), sel_tag)
          cfgs[[tag]] <- list(
            label = tag,
            desc = sprintf("PR: %d cycles, %.0f%% drop, %s",
                           cyc, drop * 100, sel_tag),
            pr_cycles = cyc,
            pr_drop = drop,
            pr_selection = sel
          )
        }
      }
    }

    # Ratchet-replacement: prune-reinsert WITH reduced ratchet
    for (cyc in c(3L, 5L)) {
      tag <- sprintf("pr_c%d_d10_noratch", cyc)
      cfgs[[tag]] <- list(
        label = tag,
        desc = sprintf("PR: %d cycles, 10%% drop, ratchet halved", cyc),
        pr_cycles = cyc,
        pr_drop = 0.10,
        pr_selection = 0L,
        ratchet_override = TRUE  # signal to halve ratchet cycles
      )
    }
  }

  cfgs
}

configs <- build_configs(stage)
total_runs <- length(configs) * length(datasets) * length(seeds)
cat(sprintf("Configs: %d, Datasets: %d, Seeds: %d -> %d total runs\n\n",
            length(configs), length(datasets), length(seeds), total_runs))

# ---- Run experiments ----
results <- data.frame(
  dataset = character(), n_tips = integer(), n_patterns = integer(),
  config = character(), seed = integer(), timeout_s = integer(),
  score = numeric(), n_trees = integer(), replicates = integer(),
  hits = integer(), wall_s = numeric(),
  pr_cycles = integer(), pr_drop = numeric(), pr_selection = integer(),
  ratchet_halved = logical(),
  tnt_best = numeric(), gap = numeric(),
  stringsAsFactors = FALSE
)

run_idx <- 0L
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  cat(sprintf("\n--- Config: %s (%s) ---\n", cfg$label, cfg$desc))

  for (ds_name in bench_names) {
    ds <- datasets[[ds_name]]
    ntip <- length(ds)
    npat <- sum(attr(ds, "weight"))

    for (s in seeds) {
      run_idx <- run_idx + 1L
      cat(sprintf("  [%d/%d] %s / %s / seed=%d ... ",
                  run_idx, total_runs, ds_name, cfg$label, s))

      set.seed(s)
      t0 <- proc.time()

      tryCatch({
        if (cfg_name == "baseline") {
          res <- MaximizeParsimony(
            ds,
            maxSeconds = timeout_s,
            strategy = "auto",
            consensusStableReps = 0L,
            nniPerturbCycles = 0L,
            driftCycles = 0L,
            verbosity = 0L,
            nThreads = 1L
          )
        } else {
          # Pass prune-reinsert params as ... args so the auto preset still
          # governs everything else (ratchetCycles, wagnerStarts, etc.).
          # Using control = SearchControl(...) marks ALL fields as explicit,
          # which discards the preset — see MaximizeParsimony.R lines 532-543.
          extra_args <- list(
            ds,
            maxSeconds = timeout_s,
            strategy = "auto",
            pruneReinsertCycles = cfg$pr_cycles,
            pruneReinsertDrop = cfg$pr_drop,
            pruneReinsertSelection = cfg$pr_selection,
            consensusStableReps = 0L,
            nniPerturbCycles = 0L,
            driftCycles = 0L,
            verbosity = 0L,
            nThreads = 1L
          )

          # Ratchet-replacement mode: halve ratchet cycles
          if (isTRUE(cfg$ratchet_override)) {
            extra_args$ratchetCycles <- 6L  # halved from preset default 12
          }

          res <- do.call(MaximizeParsimony, extra_args)
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
          config = cfg$label, seed = s, timeout_s = timeout_s,
          score = best_score, n_trees = n_trees, replicates = reps,
          hits = hits, wall_s = elapsed,
          pr_cycles = if (is.null(cfg$pr_cycles)) 0L else cfg$pr_cycles,
          pr_drop = if (is.null(cfg$pr_drop)) 0 else cfg$pr_drop,
          pr_selection = if (is.null(cfg$pr_selection)) NA_integer_
                         else cfg$pr_selection,
          ratchet_halved = isTRUE(cfg$ratchet_override),
          tnt_best = tnt_ref, gap = gap,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      })
    }
  }
}

# ---- Save results ----
outfile <- file.path(
  output_dir,
  sprintf("t289_stage%d_%ds.csv", stage, timeout_s)
)
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\n=== Results written to %s (%d rows) ===\n",
            outfile, nrow(results)))

# ---- Quick summary ----
cat("\n--- Median scores by config × dataset ---\n")
agg <- aggregate(score ~ config + dataset, data = results, FUN = median)
agg_wide <- reshape(agg, direction = "wide", idvar = "config",
                    timevar = "dataset", v.names = "score")
print(agg_wide, row.names = FALSE)

cat(sprintf("\nCompleted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
