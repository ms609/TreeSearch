# bench_pr_stage5_nni.R
#
# T-289f: Prune-reinsert Stage 5 — NNI full-tree polish cost reduction
#
# DESIGNED FOR HAMILTON HPC. Do not run locally.
#
# Stage 4 conclusion: PR (TBR full polish) is disqualified for the large preset
# at 60s budget because per-cycle cost is too high (~60s at 206 tips, leaving
# 0 replicates). Stage 5 asks whether NNI full-tree polish (pruneReinsertNni=TRUE,
# ~5x cheaper at large n) restores PR's value.
#
# Hypothesis: PR's benefit comes from topological displacement, not from the
# quality of post-reinsert local search. NNI reaches a local optimum sufficient
# to identify improvements; outer-loop TBR then polishes to full convergence.
#
# Three configs:
#   baseline:  large preset, pruneReinsertCycles=0         (no PR)
#   pr_nni:    large preset, c=5, d=5%, MISSING, NNI=TRUE  (new cheap option)
#   pr_tbr:    large preset, c=5, d=5%, MISSING, NNI=FALSE (Stage 4 reference)
#
# Same 5 datasets as Stage 4 (131-206 tips, training-split MorphoBank).
#
# Grid: 5 datasets x 3 configs x 2 budgets x 10 seeds = 300 runs
# Expected wall time: ~4-6h (pr_nni ~5x faster than pr_tbr).
#
# Usage:
#   Rscript bench_pr_stage5_nni.R [output_dir]
#   output_dir: where to write CSV. Default: "."
#
# Output: t289f_pr_nni_polish.csv

suppressPackageStartupMessages({
  library(TreeSearch)
  library(TreeTools)
})

args       <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "."

cat("=== T-289f: Prune-Reinsert Stage 5 — NNI Polish ===\n")
cat(sprintf("TreeSearch %s\n", packageVersion("TreeSearch")))
cat(sprintf("Output: %s\n", output_dir))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

# ---- Dataset definitions ----
neotrans_dir <- Sys.glob("/nobackup/*/neotrans/inst/matrices")
if (length(neotrans_dir) == 0) {
  neotrans_dir <- file.path(dirname(dirname(dirname(getwd()))),
                             "neotrans", "inst", "matrices")
}
neotrans_dir <- neotrans_dir[1]
if (!dir.exists(neotrans_dir)) stop("neotrans matrices directory not found: ", neotrans_dir)

mbank_path <- Sys.glob("/nobackup/*/TreeSearch-a/dev/benchmarks/mbank_X30754.nex")
if (length(mbank_path) == 0) {
  mbank_path <- file.path(dirname(dirname(dirname(getwd()))),
                           "TreeSearch-a", "dev", "benchmarks", "mbank_X30754.nex")
}
mbank_path <- mbank_path[1]
if (!file.exists(mbank_path)) stop("mbank_X30754.nex not found")

dataset_defs <- list(
  list(key = "mbank_X30754", path = mbank_path),
  list(key = "project4133",  path = file.path(neotrans_dir, "project4133.nex")),
  list(key = "project3701",  path = file.path(neotrans_dir, "project3701.nex")),
  list(key = "project804",   path = file.path(neotrans_dir, "project804.nex")),
  list(key = "syab07205",    path = file.path(neotrans_dir, "syab07205.nex"))
)

# ---- Config grid ----
sc_baseline <- SearchControl(
  pruneReinsertCycles = 0L
)
sc_pr_nni <- SearchControl(
  pruneReinsertCycles    = 5L,
  pruneReinsertDrop      = 0.05,
  pruneReinsertSelection = 2L,   # MISSING
  pruneReinsertNni       = TRUE  # new: NNI polish instead of TBR
)
sc_pr_tbr <- SearchControl(
  pruneReinsertCycles    = 5L,
  pruneReinsertDrop      = 0.05,
  pruneReinsertSelection = 2L,   # MISSING
  pruneReinsertNni       = FALSE # Stage 4 reference: TBR full convergence
)

configs <- list(
  baseline = sc_baseline,
  pr_nni   = sc_pr_nni,
  pr_tbr   = sc_pr_tbr
)

budgets <- c(60L, 120L)
seeds   <- 1:10

# ---- Output ----
out_file <- file.path(output_dir, "t289f_pr_nni_polish.csv")
out_cols <- c("dataset", "n_tips", "n_patterns", "config", "seed", "timeout_s",
              "score", "n_trees", "replicates", "hits", "wall_s",
              "pr_cycles", "pr_nni")
write(paste(shQuote(out_cols), collapse = ","), out_file)

total_runs <- length(dataset_defs) * length(configs) * length(budgets) * length(seeds)
cat(sprintf("Total runs: %d\n\n", total_runs))
run_i <- 0L

for (ddef in dataset_defs) {
  cat(sprintf("--- Loading: %s ---\n", ddef$key))
  ds <- tryCatch(ReadAsPhyDat(ddef$path), error = function(e) {
    cat(sprintf("  ERROR loading %s: %s\n", ddef$key, e$message))
    NULL
  })
  if (is.null(ds)) next
  n_tips     <- length(ds)
  n_patterns <- sum(attr(ds, "weight"))
  cat(sprintf("  %d taxa, %d patterns\n\n", n_tips, n_patterns))

  for (budget in budgets) {
    for (cfg_name in names(configs)) {
      sc <- configs[[cfg_name]]
      for (seed in seeds) {
        run_i <- run_i + 1L
        cat(sprintf("[%d/%d] %s | %s | budget=%ds | seed=%d ... ",
                    run_i, total_runs, ddef$key, cfg_name, budget, seed))
        t0 <- proc.time()[["elapsed"]]

        res <- tryCatch(
          MaximizeParsimony(
            dataset    = ds,
            maxSeconds = budget,
            nThreads   = 2L,
            seed       = seed,
            verbosity  = 0L,
            control    = sc
          ),
          error = function(e) {
            cat(sprintf("ERROR: %s\n", e$message))
            NULL
          }
        )

        wall_s <- proc.time()[["elapsed"]] - t0
        if (is.null(res)) next

        score      <- attr(res, "score")
        n_trees    <- length(res)
        replicates <- attr(res, "replicates")
        hits       <- attr(res, "hits")
        pr_cycles  <- if (!is.null(sc$pruneReinsertCycles)) sc$pruneReinsertCycles else 0L
        pr_nni_val <- if (!is.null(sc$pruneReinsertNni))    as.integer(sc$pruneReinsertNni) else 0L

        row <- sprintf('%s,%d,%d,%s,%d,%d,%g,%d,%d,%d,%.3f,%d,%d',
                        shQuote(ddef$key), n_tips, n_patterns, shQuote(cfg_name),
                        seed, budget, score, n_trees, replicates, hits, wall_s,
                        pr_cycles, pr_nni_val)
        write(row, out_file, append = TRUE)
        cat(sprintf("score=%g  reps=%d  wall=%.1fs\n", score, replicates, wall_s))
      }
    }
  }
}

cat(sprintf("\nDone. %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
