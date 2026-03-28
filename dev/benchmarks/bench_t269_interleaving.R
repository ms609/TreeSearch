#!/usr/bin/env Rscript
# T-269: Fine-grained sectorial interleaving benchmark
#
# DESIGNED FOR HAMILTON HPC. Do not run locally (hours of wall time).
#
# Tests whether fine-grained interleaving of sectorial search with ratchet
# perturbation improves score quality. The key question: does performing
# one sectorial pass per ratchet cycle (outerCycles = ratchetCycles) help
# compared to the current thorough preset (outerCycles = 2)?
#
# Design:
#   - Thorough preset as base (ratchetCycles=20, XSS+RSS+CSS, outerCycles=2)
#   - Vary outerCycles ∈ {1, 2, 4, 10, 20} while holding ratchetCycles=20
#   - 4 standard gap datasets (37–88 tips), 5 seeds, 30s + 60s budgets
#   - EW scoring throughout (inapplicable → missing via fitch_mode)
#
# outerCycles=1:  all 20 ratchet cycles in one block, then 1 sectorial pass
# outerCycles=2:  2 × 10 ratchet + 2 sectorial passes (current thorough)
# outerCycles=4:  4 × 5  ratchet + 4 sectorial passes
# outerCycles=10: 10 × 2 ratchet + 10 sectorial passes
# outerCycles=20: 20 × 1 ratchet + 20 sectorial passes (TNT pattern)
#
# Usage:
#   Rscript bench_t269_interleaving.R [timeout_s] [output_dir]
#   timeout_s:  search budget in seconds. Default: 30
#   output_dir: where to write CSV results. Default: "."
#
# Output: t269_interleaving_{timeout}s.csv

library(TreeSearch)
library(TreeTools)

args <- commandArgs(trailingOnly = TRUE)
timeout_s  <- if (length(args) >= 1) as.integer(args[1]) else 30L
output_dir <- if (length(args) >= 2) args[2] else "."

cat("=== T-269: Fine-Grained Sectorial Interleaving Benchmark ===\n")
cat(sprintf("Timeout: %ds\n", timeout_s))
cat(sprintf("TreeSearch version: %s\n", packageVersion("TreeSearch")))
cat(sprintf("Output dir: %s\n", output_dir))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

# ---- Datasets ----
# 4 standard datasets with persistent TNT gaps — range 37–88 tips.
# inapplicable converted to missing for EW Fitch (match TNT).
fitch_mode <- function(dataset) {
  contrast  <- attr(dataset, "contrast")
  levels    <- attr(dataset, "levels")
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

bench_names <- c("Wortley2006", "Agnarsson2004", "Zhu2013", "Dikow2009")
datasets <- lapply(
  setNames(bench_names, bench_names),
  function(nm) fitch_mode(inapplicable.phyData[[nm]])
)

# TNT reference scores (EW Fitch mode, from T-265)
tnt_best <- c(
  Wortley2006 = 479, Agnarsson2004 = 718,
  Zhu2013 = 624,     Dikow2009 = 1603
)

seeds <- 1:5

# ---- Configs ----
# Fixed thorough-preset parameters (ratchetCycles=20, no drift, no NNI-perturb)
# outerCycles varies: 1, 2, 4, 10, 20.
outer_cycles_grid <- c(1L, 2L, 4L, 10L, 20L)

build_control <- function(outer_cycles) {
  SearchControl(
    # Thorough preset base
    ratchetCycles          = 20L,
    ratchetPerturbProb     = 0.25,
    ratchetPerturbMode     = 2L,
    ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive        = FALSE,   # off for cleaner comparison
    # Vary this:
    outerCycles            = outer_cycles,
    # Sectorial
    xssRounds              = 5L,
    rssRounds              = 5L,
    cssRounds              = 2L,
    # No drift/NNI-perturb
    driftCycles            = 0L,
    nniPerturbCycles       = 0L,
    # Other thorough settings
    wagnerStarts           = 3L,
    nniFirst               = TRUE,
    consensusStableReps    = 0L
  )
}

configs <- setNames(
  lapply(outer_cycles_grid, build_control),
  sprintf("outer_%02d", outer_cycles_grid)
)

total_runs <- length(configs) * length(datasets) * length(seeds)
cat(sprintf("Configs: %d (outerCycles: %s), Datasets: %d, Seeds: %d -> %d total runs\n\n",
            length(configs),
            paste(outer_cycles_grid, collapse = "/"),
            length(datasets), length(seeds), total_runs))

# ---- Run experiments ----
results <- data.frame(
  dataset = character(), n_tips = integer(), n_patterns = integer(),
  outer_cycles = integer(), seed = integer(), timeout_s = integer(),
  score = numeric(), n_trees = integer(), replicates = integer(),
  wall_s = numeric(), tnt_best = numeric(), gap = numeric(),
  stringsAsFactors = FALSE
)

run_idx <- 0L
for (cfg_name in names(configs)) {
  ctrl <- configs[[cfg_name]]
  oc <- ctrl$outerCycles
  cat(sprintf("\n--- outerCycles = %d ---\n", oc))

  for (ds_name in bench_names) {
    ds    <- datasets[[ds_name]]
    ntip  <- length(ds)
    npat  <- sum(attr(ds, "weight"))

    for (s in seeds) {
      run_idx <- run_idx + 1L
      cat(sprintf("  [%d/%d] %s / oc=%d / seed=%d ... ",
                  run_idx, total_runs, ds_name, oc, s))

      set.seed(s)
      t0 <- proc.time()

      tryCatch({
        res <- MaximizeParsimony(
          ds,
          maxSeconds = timeout_s,
          control    = ctrl,
          verbosity  = 0L,
          nThreads   = 1L
        )

        elapsed    <- (proc.time() - t0)[3]
        best_score <- attr(res, "score")
        n_trees    <- length(res)
        reps       <- attr(res, "replicates")
        tnt_ref    <- tnt_best[ds_name]
        gap        <- if (!is.na(tnt_ref)) best_score - tnt_ref else NA_real_

        cat(sprintf("score=%g, gap=%s, reps=%d, %.1fs\n",
                    best_score,
                    if (is.na(gap)) "?" else sprintf("%+d", gap),
                    reps, elapsed))

        results <- rbind(results, data.frame(
          dataset     = ds_name, n_tips = ntip, n_patterns = npat,
          outer_cycles = oc,      seed = s,     timeout_s = timeout_s,
          score       = best_score, n_trees = n_trees, replicates = reps,
          wall_s      = elapsed,
          tnt_best    = tnt_ref,  gap = gap,
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
  sprintf("t269_interleaving_%ds.csv", timeout_s)
)
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\n=== Results written to %s (%d rows) ===\n",
            outfile, nrow(results)))

# ---- Quick summary ----
cat("\n--- Median gap by outerCycles × dataset ---\n")
agg <- aggregate(gap ~ outer_cycles + dataset, data = results, FUN = median,
                 na.rm = TRUE)
agg_wide <- reshape(agg, direction = "wide", idvar = "outer_cycles",
                    timevar = "dataset", v.names = "gap")
names(agg_wide) <- sub("gap\\.", "", names(agg_wide))
print(agg_wide[order(agg_wide$outer_cycles), ], row.names = FALSE)

cat("\n--- Median gap by outerCycles (pooled) ---\n")
agg2 <- aggregate(gap ~ outer_cycles, data = results, FUN = median,
                  na.rm = TRUE)
print(agg2[order(agg2$outer_cycles), ], row.names = FALSE)

cat(sprintf("\nCompleted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
