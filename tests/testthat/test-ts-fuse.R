library("TreeTools")

# Access unexported functions from the TreeSearch namespace
ts_fitch_score <- TreeSearch:::ts_fitch_score
ts_tbr_search <- TreeSearch:::ts_tbr_search
ts_tree_fuse <- TreeSearch:::ts_tree_fuse

# Helper: prepare dataset for ts_* functions from a phyDat object
make_ts_data <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  weight <- at$weight
  levels <- at$levels
  list(contrast = contrast, tip_data = tip_data,
       weight = weight, levels = levels)
}

# Helper: score a tree with ts engine
ts_score <- function(tree, ds) {
  ts_fitch_score(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels)
}

# Helper: run TBR search
ts_tbr <- function(tree, ds, maxHits = 1L) {
  ts_tbr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                maxHits = maxHits)
}

# Helper: run fuse
ts_fuse <- function(tree, ds, pool_edges, pool_scores,
                    accept_equal = FALSE, max_rounds = 10L) {
  ts_tree_fuse(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
               pool_edges, pool_scores,
               accept_equal = accept_equal, max_rounds = max_rounds)
}

# Helper: generate a pool of TBR-optimized trees from random starting points
make_pool <- function(ds, n_tip, n_trees, seed) {
  set.seed(seed)
  pool_edges <- vector("list", n_trees)
  pool_scores <- numeric(n_trees)
  for (i in seq_len(n_trees)) {
    start_tree <- as.phylo(sample.int(1e6, 1), n_tip)
    result <- ts_tbr(start_tree, ds, maxHits = 3L)
    pool_edges[[i]] <- result$edge
    pool_scores[i] <- result$score
  }
  list(edges = pool_edges, scores = pool_scores)
}


test_that("Fuse result has valid structure", {
  tree <- as.phylo(42, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
    1, 0, 0, 0, 1, 1, 1, 1, 0, 0
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Build a small pool
  pool <- make_pool(ds, 10, 3, seed = 8471L)

  best_idx <- which.min(pool$scores)
  result <- ts_fuse(tree, ds, pool$edges, pool$scores)

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("n_exchanges" %in% names(result))
  expect_true("n_rounds" %in% names(result))
  expect_true(result$score > 0)
})

test_that("Fusing identical trees changes nothing", {
  tree <- as.phylo(1, 12)
  set.seed(3847)
  mat <- matrix(sample(0:1, 12 * 6, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Optimize once
  tbr_result <- ts_tbr(tree, ds, maxHits = 3L)
  opt_tree <- tree
  opt_tree$edge <- tbr_result$edge
  opt_score <- ts_score(opt_tree, ds)

  # Pool of identical trees
  pool_edges <- replicate(3, tbr_result$edge, simplify = FALSE)
  pool_scores <- rep(opt_score, 3)

  result <- ts_fuse(opt_tree, ds, pool_edges, pool_scores)

  # Score should not worsen
  expect_equal(result$score, opt_score)
})

test_that("Fuse finds score <= best individual tree (20 tips)", {
  set.seed(6192)
  mat <- matrix(sample(0:1, 20 * 12, replace = TRUE),
                nrow = 20,
                dimnames = list(paste0("t", 1:20), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  pool <- make_pool(ds, 20, 5, seed = 7043L)
  best_pool_score <- min(pool$scores)
  best_idx <- which.min(pool$scores)

  start_tree <- as.phylo(1, 20)
  start_tree$edge <- pool$edges[[best_idx]]

  result <- ts_fuse(start_tree, ds, pool$edges, pool$scores)

  expect_true(result$score <= best_pool_score,
              info = paste("Fused:", result$score, "Best pool:", best_pool_score))
})

test_that("Fuse on 50-tip dataset finds at least as good as best individual", {
  skip_on_cran()
  set.seed(2518)
  mat <- matrix(sample(0:2, 50 * 20, replace = TRUE),
                nrow = 50,
                dimnames = list(paste0("t", 1:50), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  pool <- make_pool(ds, 50, 10, seed = 5839L)
  best_pool_score <- min(pool$scores)
  best_idx <- which.min(pool$scores)

  start_tree <- as.phylo(1, 50)
  start_tree$edge <- pool$edges[[best_idx]]

  result <- ts_fuse(start_tree, ds, pool$edges, pool$scores)

  expect_true(result$score <= best_pool_score,
              info = paste("Fused:", result$score, "Best pool:", best_pool_score))
})

test_that("Equal-score exchanges may change topology", {
  set.seed(4326)
  mat <- matrix(sample(0:1, 15 * 6, replace = TRUE),
                nrow = 15,
                dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  pool <- make_pool(ds, 15, 5, seed = 1958L)
  best_pool_score <- min(pool$scores)
  best_idx <- which.min(pool$scores)

  start_tree <- as.phylo(1, 15)
  start_tree$edge <- pool$edges[[best_idx]]

  result_strict <- ts_fuse(start_tree, ds, pool$edges, pool$scores,
                           accept_equal = FALSE)
  result_equal <- ts_fuse(start_tree, ds, pool$edges, pool$scores,
                          accept_equal = TRUE)

  # Both should produce valid scores
  expect_true(result_strict$score <= best_pool_score)
  expect_true(result_equal$score <= best_pool_score)
})

test_that("Fused tree score matches independent verification", {
  set.seed(9637)
  mat <- matrix(sample(0:1, 18 * 10, replace = TRUE),
                nrow = 18,
                dimnames = list(paste0("t", 1:18), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  pool <- make_pool(ds, 18, 4, seed = 3091L)
  best_idx <- which.min(pool$scores)

  start_tree <- as.phylo(1, 18)
  start_tree$edge <- pool$edges[[best_idx]]

  result <- ts_fuse(start_tree, ds, pool$edges, pool$scores)

  # Verify score independently
  result_tree <- start_tree
  result_tree$edge <- result$edge
  verified_score <- ts_score(result_tree, ds)
  expect_equal(result$score, verified_score)
})

test_that("Fuse result tree has valid topology", {
  tree <- as.phylo(100, 14)
  set.seed(2765)
  mat <- matrix(sample(0:1, 14 * 6, replace = TRUE),
                nrow = 14,
                dimnames = list(paste0("t", 1:14), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  pool <- make_pool(ds, 14, 3, seed = 8122L)
  best_idx <- which.min(pool$scores)

  start_tree <- tree
  start_tree$edge <- pool$edges[[best_idx]]

  result <- ts_fuse(start_tree, ds, pool$edges, pool$scores)

  # Right number of edges
  n_tip <- length(tree$tip.label)
  expect_equal(nrow(result$edge), 2 * (n_tip - 1))

  # All tips present
  tips_in_tree <- sort(result$edge[result$edge[, 2] <= n_tip, 2])
  expect_equal(tips_in_tree, seq_len(n_tip))
})

test_that("Fuse exchanges clades when bipartition has flipped orientation", {
  # Two 8-tip trees share bipartition {t3,t4,t5,t6}|{t1,t2,t7,t8}.
  # In tree A, the node rooting {t3,t4,t5,t6} has tip 0 (t1) outside (was_flipped=false).
  # In tree B, the matching node roots {t1,t2,t7,t8} with tip 0 inside (was_flipped=true).
  # Without re-rooting, fusing skips this bipartition entirely.
  # With re-rooting at tip 0, both sides are consistently oriented and the

  # exchange proceeds.

  # Tree A: (t1, ((t2,(t7,t8)), ((t3,t4),(t5,t6))))
  edgeA <- matrix(c(
    9, 1,   9, 10,  10, 11,  10, 13,  11, 2,  11, 12,
    12, 7,  12, 8,  13, 14,  13, 15,  14, 3,  14, 4,
    15, 5,  15, 6
  ), ncol = 2, byrow = TRUE)
  treeA <- structure(list(edge = edgeA, tip.label = paste0("t", 1:8),
                          Nnode = 7L), class = "phylo")
  treeA <- Preorder(treeA)

  # Tree B: (t3, (((t1,t2),(t7,t8)), (t4,(t5,t6))))
  edgeB <- matrix(c(
    9, 3,   9, 10,  10, 11,  10, 14,  11, 12,  11, 13,
    12, 1,  12, 2,  13, 7,   13, 8,   14, 4,   14, 15,
    15, 5,  15, 6
  ), ncol = 2, byrow = TRUE)
  treeB <- structure(list(edge = edgeB, tip.label = paste0("t", 1:8),
                          Nnode = 7L), class = "phylo")
  treeB <- Preorder(treeB)

  # Dataset where tree B's arrangement is strictly better.
  # t4,t5,t6 form a nested group (favors (t3,(t4,(t5,t6))) over ((t3,t4),(t5,t6)))
  mat <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0,  # t1
    0, 0, 0, 0, 0, 1, 0, 0,  # t2
    1, 0, 0, 0, 0, 0, 0, 0,  # t3
    1, 1, 1, 0, 0, 0, 0, 0,  # t4
    1, 1, 1, 1, 0, 0, 0, 0,  # t5
    1, 1, 1, 1, 0, 0, 0, 0,  # t6
    0, 0, 0, 0, 0, 0, 1, 0,  # t7
    0, 0, 0, 0, 0, 0, 1, 1   # t8
  ), nrow = 8, byrow = TRUE, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  scoreA <- ts_score(treeA, ds)
  scoreB <- ts_score(treeB, ds)
  expect_true(scoreB < scoreA,
              info = "Donor should score strictly better than recipient")

  result <- ts_fuse(treeA, ds, list(treeB$edge), c(scoreB))

  # With re-rooting, fusing should find the exchange and improve the score
  expect_true(result$n_exchanges > 0,
              info = "Should find exchange via re-rooted bipartition matching")
  expect_true(result$score <= scoreB,
              info = paste("Fused:", result$score, "Donor:", scoreB))
})
