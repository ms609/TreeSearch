library("TreeTools")

test_that("Pool deduplication: same tree added twice → pool size 1", {
  tree <- as.phylo(42, 8)
  edges <- list(tree$edge, tree$edge)
  scores <- c(10, 10)

  result <- ts_pool_test(edges, scores, 8L)

  expect_true(result$added[1])
  expect_false(result$added[2])
  expect_equal(result$pool_size, 1)
})

test_that("Pool stores different trees", {
  tree1 <- as.phylo(1, 8)
  tree2 <- as.phylo(2, 8)
  edges <- list(tree1$edge, tree2$edge)
  scores <- c(10, 10)

  result <- ts_pool_test(edges, scores, 8L)

  expect_true(result$added[1])
  expect_true(result$added[2])
  expect_equal(result$pool_size, 2)
})

test_that("Pool tracks best score", {
  tree1 <- as.phylo(1, 8)
  tree2 <- as.phylo(2, 8)
  edges <- list(tree1$edge, tree2$edge)
  scores <- c(15, 10)

  result <- ts_pool_test(edges, scores, 8L)

  expect_equal(result$best_score, 10)
})

test_that("Pool evicts trees beyond suboptimal threshold", {
  tree1 <- as.phylo(1, 8)
  tree2 <- as.phylo(2, 8)
  tree3 <- as.phylo(3, 8)
  edges <- list(tree1$edge, tree2$edge, tree3$edge)
  # tree1 score 15, tree2 score 10, tree3 score 10
  # With suboptimal = 2, tree1 (15 > 10+2 = 12) should be evicted
  scores <- c(15, 10, 10)

  result <- ts_pool_test(edges, scores, 8L, max_size = 100L,
                         suboptimal = 2.0)

  # tree1 added first (best so far), tree2 triggers eviction of tree1

  expect_equal(result$best_score, 10)
  expect_equal(result$pool_size, 2)  # tree2 and tree3 survive
})

test_that("hits_to_best counts rediscoveries", {
  tree1 <- as.phylo(1, 8)
  tree2 <- as.phylo(2, 8)
  edges <- list(tree1$edge, tree2$edge, tree1$edge)
  scores <- c(10, 10, 10)

  result <- ts_pool_test(edges, scores, 8L)

  # tree1: added, hits=1
  # tree2: added (different topology, same score), hits=2
  # tree1 again: duplicate, hits=3
  expect_equal(result$hits_to_best, 3)
})

test_that("Pool respects max_size", {
  n_trees <- 5
  trees <- lapply(seq_len(n_trees), function(i) as.phylo(i, 8))
  edges <- lapply(trees, `[[`, "edge")
  scores <- rep(10, n_trees)

  result <- ts_pool_test(edges, scores, 8L, max_size = 3L)

  expect_true(result$pool_size <= 3)
})

test_that("Pool rejects trees worse than threshold", {
  tree1 <- as.phylo(1, 8)
  tree2 <- as.phylo(2, 8)
  edges <- list(tree1$edge, tree2$edge)
  # tree1 at score 10, tree2 at score 20 with suboptimal=5
  # tree2 (20 > 10+5=15) should be rejected
  scores <- c(10, 20)

  result <- ts_pool_test(edges, scores, 8L, max_size = 100L,
                         suboptimal = 5.0)

  expect_equal(result$pool_size, 1)
  expect_true(result$added[1])
  expect_false(result$added[2])
})

test_that("Deep copy: modifying original doesn't corrupt pool", {
  # This tests that TreeState copies in the pool are independent.
  # We test indirectly: add a tree, then add a different tree,

  # and verify the pool still has correct topology for both.
  tree1 <- as.phylo(1, 8)
  tree2 <- as.phylo(2, 8)

  # If the pool stored shallow copies, both entries would point to the same

  # data, and adding tree2 would overwrite tree1's topology.
  result <- ts_pool_test(list(tree1$edge, tree2$edge), c(10, 10), 8L)
  expect_equal(result$pool_size, 2)

  # Verify they're actually different by checking ts_trees_equal
  expect_false(ts_trees_equal(tree1$edge, tree2$edge, 8))
})
