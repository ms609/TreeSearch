#!/usr/bin/env Rscript
# T-289d: Prune-reinsert Stage 3 — new drop criteria (MISSING, COMBINED)
#
# DESIGNED FOR HAMILTON HPC. Do not run locally.
#
# Stage 2 (9 configs x 10 seeds x 60s, mbank_X30754) established:
#   - All PR configs improve over baseline at 180t.
#   - Instability-weighted dropping (sel=1) beats random (sel=0) by 1.8–3.3 steps.
#   - pr_c5_d05 (−12.3 steps, 3.0 reps) best cost-quality ratio at sel=0.
#   - pr_c5_d10_sel1 (−14.1 steps, 2.2 reps) best overall.
#   - Gap: pr_c5_d05_sel1 not tested.
#
# Stage 3 goals:
#   1. Fill gap: pr_c5_d05_sel1 (instability-weighted at cheapest good config).
#   2. Benchmark new criteria: MISSING (sel=2), COMBINED (sel=3) at d05 and d10.
#   3. Reference repeats: baseline + pr_c5_d05_sel0 + pr_c5_d10_sel1 for
#      within-run comparability (avoids cross-run seed variance).
#
# Grid: 8 configs × 1 dataset × 10 seeds × 60s ≈ 87 min wall time.
#
# Drop criteria (pruneReinsertSelection):
#   0 = RANDOM       uniform random
#   1 = INSTABILITY  weighted by positional instability in pool
#   2 = MISSING      weighted by ambiguous/inapplicable character count
#   3 = COMBINED     instability × (1 + normalised missingness)
#
# Usage:
#   Rscript bench_pr_stage3_mbank.R [timeout_s] [output_dir]

library(TreeSearch)
library(TreeTools)

args       <- commandArgs(trailingOnly = TRUE)
timeout_s  <- if (length(args) >= 1) as.integer(args[1]) else 60L
output_dir <- if (length(args) >= 2) args[2] else "."

cat("=== T-289d: Prune-Reinsert Stage 3 (new criteria, mbank, Brazeau) ===\n")
cat(sprintf("Timeout: %ds  |  TreeSearch %s\n", timeout_s,
            packageVersion("TreeSearch")))
cat(sprintf("Output: %s\n", output_dir))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

# ---- Load 180-tip dataset ----
mbank_path <- Sys.glob("/nobackup/*/TreeSearch-a/dev/benchmarks/mbank_X30754.nex")
if (length(mbank_path) == 0) {
  mbank_path <- file.path(dirname(dirname(dirname(getwd()))),
                          "TreeSearch-a", "dev", "benchmarks", "mbank_X30754.nex")
}
if (length(mbank_path) > 0) mbank_path <- mbank_path[1]
if (!file.exists(mbank_path)) stop("mbank_X30754.nex not found")
cat("Loading:", mbank_path, "\n")
ds <- ReadAsPhyDat(mbank_path)
cat(sprintf("  %d taxa, %d patterns\n\n", length(ds), sum(attr(ds, "weight"))))

seeds <- 1:10

# ---- Config grid ----
#
# Notation: pr_c{cycles}_d{drop%}_sel{selection}
# References from Stage 2 included for within-run comparability.
configs <- list(
  baseline = list(
    label = "baseline", desc = "No prune-reinsert",
    pr_cycles = 0L, pr_drop = 0.0, pr_selection = 0L
  ),
  # --- d=5%, c=5: cheapest good config from Stage 2 ---
  pr_c5_d05_sel0 = list(
    label = "pr_c5_d05_sel0", desc = "c5 d5% random (Stage2 ref)",
    pr_cycles = 5L, pr_drop = 0.05, pr_selection = 0L
  ),
  pr_c5_d05_sel1 = list(
    label = "pr_c5_d05_sel1", desc = "c5 d5% instability (gap)",
    pr_cycles = 5L, pr_drop = 0.05, pr_selection = 1L
  ),
  pr_c5_d05_sel2 = list(
    label = "pr_c5_d05_sel2", desc = "c5 d5% missing (new)",
    pr_cycles = 5L, pr_drop = 0.05, pr_selection = 2L
  ),
  pr_c5_d05_sel3 = list(
    label = "pr_c5_d05_sel3", desc = "c5 d5% combined (new)",
    pr_cycles = 5L, pr_drop = 0.05, pr_selection = 3L
  ),
  # --- d=10%, c=5: Stage 2 overall winner config ---
  pr_c5_d10_sel1 = list(
    label = "pr_c5_d10_sel1", desc = "c5 d10% instability (Stage2 ref)",
    pr_cycles = 5L, pr_drop = 0.10, pr_selection = 1L
  ),
  pr_c5_d10_sel2 = list(
    label = "pr_c5_d10_sel2", desc = "c5 d10% missing (new)",
    pr_cycles = 5L, pr_drop = 0.10, pr_selection = 2L
  ),
  pr_c5_d10_sel3 = list(
    label = "pr_c5_d10_sel3", desc = "c5 d10% combined (new)",
    pr_cycles = 5L, pr_drop = 0.10, pr_selection = 3L
  )
)

total_runs <- length(configs) * length(seeds)
cat(sprintf("Configs: %d, Seeds: %d -> %d total runs\n\n",
            length(configs), length(seeds), total_runs))

# ---- Run experiments ----
results <- data.frame(
  dataset = character(), n_tips = integer(), n_patterns = integer(),
  config = character(), seed = integer(), timeout_s = integer(),
  score = numeric(), n_trees = integer(), replicates = integer(),
  hits = integer(), wall_s = numeric(),
  pr_cycles = integer(), pr_drop = numeric(), pr_selection = integer(),
  stringsAsFactors = FALSE
)

ntip <- length(ds)
npat <- sum(attr(ds, "weight"))
run_idx <- 0L
outfile  <- file.path(output_dir, sprintf("t289d_stage3_%ds.csv", timeout_s))

for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  cat(sprintf("\n--- %s: %s ---\n", cfg$label, cfg$desc))

  for (s in seeds) {
    run_idx <- run_idx + 1L
    cat(sprintf("  [%d/%d] seed=%d ... ", run_idx, total_runs, s))

    set.seed(s)
    t0 <- proc.time()

    tryCatch({
      if (cfg$pr_cycles == 0L) {
        res <- MaximizeParsimony(
          ds,
          maxSeconds           = timeout_s,
          strategy             = "auto",
          consensusStableReps  = 0L,
          nniPerturbCycles     = 0L,
          driftCycles          = 0L,
          verbosity            = 0L,
          nThreads             = 1L
        )
      } else {
        res <- MaximizeParsimony(
          ds,
          maxSeconds             = timeout_s,
          strategy               = "auto",
          pruneReinsertCycles    = cfg$pr_cycles,
          pruneReinsertDrop      = cfg$pr_drop,
          pruneReinsertSelection = cfg$pr_selection,
          consensusStableReps    = 0L,
          nniPerturbCycles       = 0L,
          driftCycles            = 0L,
          verbosity              = 0L,
          nThreads               = 1L
        )
      }

      elapsed    <- (proc.time() - t0)[3]
      best_score <- attr(res, "score")
      reps       <- attr(res, "replicates")
      hits       <- attr(res, "hits")

      cat(sprintf("score=%g, reps=%d, %.1fs\n", best_score, reps, elapsed))

      results <- rbind(results, data.frame(
        dataset = "mbank_X30754", n_tips = ntip, n_patterns = npat,
        config = cfg$label, seed = s, timeout_s = timeout_s,
        score = best_score, n_trees = length(res), replicates = reps,
        hits = hits, wall_s = elapsed,
        pr_cycles = cfg$pr_cycles, pr_drop = cfg$pr_drop,
        pr_selection = cfg$pr_selection,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    })
  }

  # Save after each config (crash recovery)
  write.csv(results, outfile, row.names = FALSE)
}

write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\n=== Results written to %s (%d rows) ===\n", outfile, nrow(results)))

# ---- Quick summary ----
cat("\n--- Mean delta vs baseline ---\n")
bl_mean <- mean(results$score[results$config == "baseline"])
agg <- aggregate(score ~ config + pr_selection, data = results, FUN = mean)
agg$delta <- round(agg$score - bl_mean, 2)
agg <- agg[order(agg$delta), ]
print(agg, row.names = FALSE)

cat(sprintf("\nCompleted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
