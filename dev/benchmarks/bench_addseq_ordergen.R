# Order-generation wall-cost: how much does building N Wagner trees under
# each addition-sequence criterion cost in isolation (no downstream TBR)?
#
# This is the other half of the mission question. bench_addseq_cell.R
# measures REACHABILITY (score vs replicate index, deterministic/load-
# invariant); this script measures the fixed per-replicate WALL-CLOCK
# overhead the extra order-generation work adds (the O(n^2) distance
# matrix / per-character bookkeeping for CLOSEST/FURTHEST/INFORMATIVE vs
# RANDOM's O(n) Fisher-Yates shuffle). Net verdict = reachability gain minus
# this cost.
#
# Must run standalone / single-tenant: fast-iteration.md is explicit that
# any wall-clock reading is only trustworthy "with the pool drained" (a
# shared, contended node invalidates timing, though not the deterministic
# reachability curve from the other script). Submit this as its own small
# single-task job, not part of the array.
#
# Env: TS_LIB, TS_DATASETS, TS_BIASES, TS_NREPS, PARTIAL_DIR.
# Local test: Rscript dev/benchmarks/bench_addseq_ordergen.R

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})

n_reps <- as.integer(Sys.getenv("TS_NREPS", "50"))
biases <- as.integer(strsplit(trimws(Sys.getenv("TS_BIASES", "0 3 4 5")), "\\s+")[[1]])
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
        "Longrich2010 Vinther2008 Sansom2010 Wortley2006 Eklund2004 Zanol2014 Zhu2013 Dikow2009")),
        "\\s+")[[1]]
outdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_addseq")

data("inapplicable.phyData", package = "TreeSearch")

results <- list()
for (nm in dsN) {
  pd <- inapplicable.phyData[[nm]]
  if (is.null(pd)) { warning("Dataset not found: ", nm); next }
  at <- attributes(pd)
  contrast <- at$contrast
  tip_data <- matrix(unlist(pd, use.names = FALSE), nrow = length(pd), byrow = TRUE)
  weight <- TreeSearch:::.ScaleWeight(at$weight)
  levels <- at$levels

  for (b in biases) {
    temp <- if (b == 0L) 1.0 else 0.3
    set.seed(1)
    t0 <- proc.time()[["elapsed"]]
    res <- TreeSearch:::ts_wagner_bias_bench(
      contrast, tip_data, weight, levels, integer(0), -1.0,
      b, temp, n_reps, FALSE  # run_tbr = FALSE: time construction only
    )
    wall <- proc.time()[["elapsed"]] - t0
    results[[length(results) + 1L]] <- data.frame(
      dataset = nm, bias = b, n_reps = n_reps,
      total_wall_s = wall, per_tree_ms = 1000 * wall / n_reps,
      mean_wagner_score = mean(res$wagner_score),
      sd_wagner_score = sd(res$wagner_score),
      stringsAsFactors = FALSE
    )
    cat(sprintf("%-14s bias=%d: %7.2f ms/tree (n=%d, mean score=%.1f)\n",
                nm, b, 1000 * wall / n_reps, n_reps, mean(res$wagner_score)))
  }
}

out <- do.call(rbind, results)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(outdir, "ordergen_wallcost.csv"), row.names = FALSE)
cat(sprintf("\nWritten %d rows -> %s/ordergen_wallcost.csv\n", nrow(out), outdir))
