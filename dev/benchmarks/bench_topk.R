#!/usr/bin/env Rscript
# Benchmark: screeningTopK = 1 vs 5
# Compare search quality (MCI) and wall time with batch CID evaluation.

suppressPackageStartupMessages({
  library(TreeSearch)
  library(TreeTools)
})

cat("TreeSearch version:", as.character(packageVersion("TreeSearch")), "\n")

# --- Generate trees (same as treeSample benchmark) ---
generate_trees <- function(n_tip, n_tree, seed) {
  set.seed(seed)
  base <- RandomTree(n_tip, root = TRUE)
  n_internal <- n_tip - 1L
  trees <- vector("list", n_tree)
  trees[[1]] <- base
  for (i in 2:n_tree) {
    tr <- base
    n_perturb <- sample(1:5, 1)
    for (j in seq_len(n_perturb)) {
      int_nodes <- unique(tr$edge[, 1])
      node <- sample(int_nodes, 1)
      children <- tr$edge[tr$edge[, 1] == node, 2]
      if (length(children) == 2 && all(children > n_tip)) {
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

run_one <- function(trees, topk, max_seconds, seed) {
  set.seed(seed)
  t0 <- proc.time()
  result <- InfoConsensus(
    trees,
    treeSample = "auto",
    screeningTopK = topk,
    maxSeconds = max_seconds,
    neverDrop = TRUE,
    verbosity = 0L,
    nThreads = 1L
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  list(
    score = attr(result, "score"),
    n_tip = Ntip(result),
    elapsed = as.numeric(elapsed)
  )
}

# --- Grid ---
configs <- expand.grid(
  n_tip = c(50, 100),
  budget = c(30, 60),
  seed = c(2847, 5193, 7641),
  stringsAsFactors = FALSE
)

cat("Grid:", nrow(configs), "conditions × 2 topK levels\n")
cat("Generating trees...\n")

tree_cache <- list()
for (n_tip in unique(configs$n_tip)) {
  key <- paste0("t", n_tip)
  cat("  Generating 1000 trees with", n_tip, "tips...\n")
  tree_cache[[key]] <- generate_trees(n_tip, 1000, seed = 3917 + n_tip)
}

# --- Run ---
results <- data.frame()

for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, ]
  key <- paste0("t", cfg$n_tip)
  trees <- tree_cache[[key]]

  for (topk in c(1L, 5L)) {
    cat(sprintf("[%d/%d] n_tip=%d, budget=%ds, seed=%d, topK=%d\n",
                i, nrow(configs), cfg$n_tip, cfg$budget, cfg$seed, topk))
    tryCatch({
      res <- run_one(trees, topk, cfg$budget, cfg$seed)
      results <- rbind(results, data.frame(
        n_tip = cfg$n_tip,
        budget = cfg$budget,
        seed = cfg$seed,
        topK = topk,
        score = res$score,
        elapsed = res$elapsed,
        stringsAsFactors = FALSE
      ))
      cat(sprintf("  -> score=%.4f, elapsed=%.1fs\n", res$score, res$elapsed))
    }, error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      results <<- rbind(results, data.frame(
        n_tip = cfg$n_tip, budget = cfg$budget, seed = cfg$seed,
        topK = topk, score = NA_real_, elapsed = NA_real_,
        stringsAsFactors = FALSE
      ))
    })
  }
}

# --- Save ---
outfile <- sprintf("topk_bench_%s.csv", format(Sys.time(), "%Y%m%d_%H%M"))
write.csv(results, outfile, row.names = FALSE)
cat("\nResults saved to:", outfile, "\n")

# --- Summary ---
cat("\n=== Summary ===\n")
if (nrow(results) > 0 && !all(is.na(results$score))) {
  agg <- aggregate(
    cbind(score, elapsed) ~ n_tip + budget + topK,
    data = results,
    FUN = function(x) c(median = median(x), mean = mean(x))
  )
  print(agg)
}
