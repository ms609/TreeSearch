#!/usr/bin/env Rscript
# Benchmark: treeSample = "auto" vs treeSample = Inf
# Run on Hamilton HPC (SBATCH wrapper below)
#
# Design:
#   For each (n_tip, n_tree, budget) combination:
#   - Run InfoConsensus with treeSample = "auto" (subsampled Phase 1)
#   - Run InfoConsensus with treeSample = Inf   (full tree set)
#   Compare: wall time, MCI score, number of search replicates completed
#
# Trees: generated via bootstrap-like NNI perturbation of a random base tree,
# so input trees share partial structure (realistic split overlap).

suppressPackageStartupMessages({
  library(TreeSearch)
  library(TreeTools)
})

cat("TreeSearch version:", as.character(packageVersion("TreeSearch")), "\n")

# --- Helpers ---
generate_trees <- function(n_tip, n_tree, seed) {
  set.seed(seed)
  base <- RandomTree(n_tip, root = TRUE)
  edges <- base$edge
  n_internal <- n_tip - 1L

  trees <- vector("list", n_tree)
  trees[[1]] <- base
  for (i in 2:n_tree) {
    tr <- base
    # 1-5 random NNI-like perturbations (swap children of random internal nodes)
    n_perturb <- sample(1:5, 1)
    for (j in seq_len(n_perturb)) {
      # Pick a random internal node with 2 internal children
      int_nodes <- unique(tr$edge[, 1])
      node <- sample(int_nodes, 1)
      children <- tr$edge[tr$edge[, 1] == node, 2]
      if (length(children) == 2 && all(children > n_tip)) {
        # Swap one grandchild between children
        gc1 <- tr$edge[tr$edge[, 1] == children[1], 2]
        gc2 <- tr$edge[tr$edge[, 1] == children[2], 2]
        if (length(gc1) >= 1 && length(gc2) >= 1) {
          swap1 <- sample(gc1, 1)
          swap2 <- sample(gc2, 1)
          tr$edge[tr$edge[, 2] == swap1, 1] <- children[2]
          tr$edge[tr$edge[, 2] == swap2, 1] <- children[1]
        }
      }
    }
    tr <- Preorder(tr)
    trees[[i]] <- tr
  }
  class(trees) <- "multiPhylo"
  trees
}

run_one <- function(trees, tree_sample, max_seconds, seed) {
  set.seed(seed)
  t0 <- proc.time()
  result <- InfoConsensus(
    trees,
    treeSample = tree_sample,
    maxSeconds = max_seconds,
    neverDrop = TRUE,
    verbosity = 0L,
    nThreads = 1L
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  score <- attr(result, "score")
  n_tip_result <- Ntip(result)

  list(
    score = score,
    n_tip = n_tip_result,
    elapsed = as.numeric(elapsed),
    seed = seed
  )
}

# --- Experimental grid ---
configs <- expand.grid(
  n_tip = c(50, 100),
  n_tree = 1000,
  budget = c(30, 60, 120),
  seed = c(2847, 5193, 7641),
  stringsAsFactors = FALSE
)

cat("Experimental grid:", nrow(configs), "conditions × 2 treeSample levels\n")
cat("Generating trees...\n")

# Pre-generate trees (same trees for both auto and Inf)
tree_cache <- list()
for (n_tip in unique(configs$n_tip)) {
  key <- paste0("t", n_tip)
  cat("  Generating", 1000, "trees with", n_tip, "tips...\n")
  tree_cache[[key]] <- generate_trees(n_tip, 1000, seed = 3917 + n_tip)
  cat("    Done. Auto treeSample =", min(1000, max(50, 2 * n_tip)), "\n")
}

# --- Run benchmarks ---
results <- data.frame()

for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, ]
  key <- paste0("t", cfg$n_tip)
  trees <- tree_cache[[key]]

  for (ts in c("auto", Inf)) {
    ts_label <- if (is.character(ts)) ts else "inf"
    cat(sprintf("[%d/%d] n_tip=%d, budget=%ds, seed=%d, treeSample=%s\n",
                i, nrow(configs), cfg$n_tip, cfg$budget, cfg$seed, ts_label))

    tryCatch({
      res <- run_one(trees, ts, cfg$budget, cfg$seed)
      results <- rbind(results, data.frame(
        n_tip = cfg$n_tip,
        n_tree = cfg$n_tree,
        budget = cfg$budget,
        seed = cfg$seed,
        treeSample = ts_label,
        score = res$score,
        n_tip_result = res$n_tip,
        elapsed = res$elapsed,
        stringsAsFactors = FALSE
      ))
      cat(sprintf("  -> score=%.4f, elapsed=%.1fs\n", res$score, res$elapsed))
    }, error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      results <<- rbind(results, data.frame(
        n_tip = cfg$n_tip,
        n_tree = cfg$n_tree,
        budget = cfg$budget,
        seed = cfg$seed,
        treeSample = ts_label,
        score = NA_real_,
        n_tip_result = NA_integer_,
        elapsed = NA_real_,
        stringsAsFactors = FALSE
      ))
    })
  }
}

# --- Save results ---
outfile <- sprintf("treesample_bench_%s.csv",
                   format(Sys.time(), "%Y%m%d_%H%M"))
write.csv(results, outfile, row.names = FALSE)
cat("\nResults saved to:", outfile, "\n")

# --- Summary ---
cat("\n=== Summary ===\n")
if (nrow(results) > 0 && !all(is.na(results$score))) {
  library(stats)
  agg <- aggregate(
    cbind(score, elapsed) ~ n_tip + budget + treeSample,
    data = results,
    FUN = function(x) c(median = median(x), mean = mean(x))
  )
  print(agg)
}
