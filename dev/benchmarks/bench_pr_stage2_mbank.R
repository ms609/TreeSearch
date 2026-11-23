#!/usr/bin/env Rscript
# T-289c: Prune-reinsert Stage 2 — mbank_X30754 (180t) only, Brazeau scoring
#
# DESIGNED FOR HAMILTON HPC. Do not run locally.
#
# Stage 1 (13 configs × 5 datasets × 5 seeds × 30s) showed:
#   - ≤88t: PR is net-negative (replicate cost >> score gain). No further testing.
#   - 180t: Real signal. Best configs by mean delta vs baseline:
#       pr_c3_d10: −8.0 (4/5 seeds), pr_c5_d10: −6.6 (5/5 seeds, most consistent)
#       pr_c5_d05: −6.8 (4/5),       pr_c3_d05: −4.8 (3/5)
#       pr_c1_d10: −2.8 (3/5) — weak but cheap
#     d≥20% with c≥3 rarely completes a replicate in 30s.
#
# Stage 2 goals:
#   1. Confirm signal at 60s (≥2 completed replicates per seed).
#   2. Narrow to best cycle/drop combination.
#   3. Test selection=1 (greedy insertion) for top-2 configs.
#
# Configs tested (8 + baseline = 9 total):
#   baseline, pr_c1_d10,
#   pr_c3_d05, pr_c3_d10, pr_c3_d10_sel1,
#   pr_c5_d05, pr_c5_d10, pr_c5_d10_sel1
#
# Grid: 9 configs × 1 dataset × 10 seeds × 60s ≈ 90 min wall time.
#
# Usage:
#   Rscript bench_pr_stage2_mbank.R [timeout_s] [output_dir]
#   timeout_s:  search budget in seconds. Default: 60
#   output_dir: where to write CSV. Default: "."
#
# Output: t289c_stage2_{timeout}s.csv

library(TreeSearch)
library(TreeTools)

args       <- commandArgs(trailingOnly = TRUE)
timeout_s  <- if (length(args) >= 1) as.integer(args[1]) else 60L
output_dir <- if (length(args) >= 2) args[2] else "."

cat("=== T-289c: Prune-Reinsert Stage 2 (mbank, Brazeau) ===\n")
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
# Stage 1 top performers (all random selection, pr_selection=0):
#   pr_c3_d10:  mean delta −8.0, 4/5 seeds improved
#   pr_c5_d10:  mean delta −6.6, 5/5 seeds improved  ← most consistent
#   pr_c5_d05:  mean delta −6.8, 4/5
#   pr_c3_d05:  mean delta −4.8, 3/5
#   pr_c1_d10:  mean delta −2.8, 3/5  — cheap reference
#
# Also test selection=1 (greedy insertion) for the top-2 configs.
configs <- list(
  baseline = list(
    label = "baseline",
    desc  = "No prune-reinsert (auto preset)",
    pr_cycles = 0L, pr_drop = 0.0, pr_selection = 0L
  ),
  pr_c1_d10 = list(
    label = "pr_c1_d10",
    desc  = "PR 1 cycle, 10% drop, random",
    pr_cycles = 1L, pr_drop = 0.10, pr_selection = 0L
  ),
  pr_c3_d05 = list(
    label = "pr_c3_d05",
    desc  = "PR 3 cycles, 5% drop, random",
    pr_cycles = 3L, pr_drop = 0.05, pr_selection = 0L
  ),
  pr_c3_d10 = list(
    label = "pr_c3_d10",
    desc  = "PR 3 cycles, 10% drop, random",
    pr_cycles = 3L, pr_drop = 0.10, pr_selection = 0L
  ),
  pr_c3_d10_sel1 = list(
    label = "pr_c3_d10_sel1",
    desc  = "PR 3 cycles, 10% drop, greedy insertion",
    pr_cycles = 3L, pr_drop = 0.10, pr_selection = 1L
  ),
  pr_c5_d05 = list(
    label = "pr_c5_d05",
    desc  = "PR 5 cycles, 5% drop, random",
    pr_cycles = 5L, pr_drop = 0.05, pr_selection = 0L
  ),
  pr_c5_d10 = list(
    label = "pr_c5_d10",
    desc  = "PR 5 cycles, 10% drop, random",
    pr_cycles = 5L, pr_drop = 0.10, pr_selection = 0L
  ),
  pr_c5_d10_sel1 = list(
    label = "pr_c5_d10_sel1",
    desc  = "PR 5 cycles, 10% drop, greedy insertion",
    pr_cycles = 5L, pr_drop = 0.10, pr_selection = 1L
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
          maxSeconds      = timeout_s,
          strategy        = "auto",
          consensusStableReps = 0L,
          nniPerturbCycles = 0L,
          driftCycles     = 0L,
          verbosity       = 0L,
          nThreads        = 1L
        )
      } else {
        res <- MaximizeParsimony(
          ds,
          maxSeconds           = timeout_s,
          strategy             = "auto",
          pruneReinsertCycles  = cfg$pr_cycles,
          pruneReinsertDrop    = cfg$pr_drop,
          pruneReinsertSelection = cfg$pr_selection,
          consensusStableReps  = 0L,
          nniPerturbCycles     = 0L,
          driftCycles          = 0L,
          verbosity            = 0L,
          nThreads             = 1L
        )
      }

      elapsed <- (proc.time() - t0)[3]
      best_score <- attr(res, "score")
      n_trees    <- length(res)
      reps       <- attr(res, "replicates")
      hits       <- attr(res, "hits")

      cat(sprintf("score=%g, reps=%d, %.1fs\n", best_score, reps, elapsed))

      results <- rbind(results, data.frame(
        dataset = "mbank_X30754", n_tips = ntip, n_patterns = npat,
        config = cfg$label, seed = s, timeout_s = timeout_s,
        score = best_score, n_trees = n_trees, replicates = reps,
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
  outfile <- file.path(output_dir,
                       sprintf("t289c_stage2_%ds.csv", timeout_s))
  write.csv(results, outfile, row.names = FALSE)
}

# ---- Save final ----
outfile <- file.path(output_dir, sprintf("t289c_stage2_%ds.csv", timeout_s))
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\n=== Results written to %s (%d rows) ===\n",
            outfile, nrow(results)))

# ---- Quick summary ----
cat("\n--- Mean scores by config ---\n")
agg <- aggregate(score ~ config, data = results, FUN = mean)
bl  <- agg$score[agg$config == "baseline"]
agg$delta <- round(agg$score - bl, 2)
agg <- agg[order(agg$delta), ]
print(agg, row.names = FALSE)

cat(sprintf("\nCompleted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
