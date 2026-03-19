# Tier 1: runs on CRAN
# Tests for WideSample() — greedy maximin tree subsampling

test_that("n >= length(trees) returns all trees", {
  trees <- as.phylo(0:9, nTip = 8)
  attr(trees, "score") <- 42
  result <- WideSample(trees, 10)
  expect_equal(length(result), 10)
  expect_equal(attr(result, "score"), 42)
  result2 <- WideSample(trees, 20)
  expect_equal(length(result2), 10)
})

test_that("n = 0 returns empty multiPhylo", {
  trees <- as.phylo(0:4, nTip = 8)
  result <- WideSample(trees, 0)
  expect_s3_class(result, "multiPhylo")
  expect_equal(length(result), 0)
})

test_that("n = 1 returns a single tree", {
  trees <- as.phylo(0:9, nTip = 8)
  set.seed(4821)
  result <- WideSample(trees, 1)
  expect_equal(length(result), 1)
})

test_that("random method returns correct size", {
  trees <- as.phylo(0:49, nTip = 8)
  set.seed(7193)
  result <- WideSample(trees, 10, method = "random")
  expect_equal(length(result), 10)
})

test_that("maximin selects more diverse trees than random", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:99, nTip = 12)
  n <- 10

  set.seed(5028)
  wide <- WideSample(trees, n, method = "maximin")
  # Compare minimum pairwise distance across many random samples
  min_dist_wide <- min(as.matrix(
    TreeDist::ClusteringInfoDistance(wide)
  )[lower.tri(diag(n))])

  set.seed(5028)
  min_dists_random <- vapply(1:20, function(i) {
    r <- WideSample(trees, n, method = "random")
    min(as.matrix(
      TreeDist::ClusteringInfoDistance(r)
    )[lower.tri(diag(n))])
  }, double(1))

  # Maximin's minimum gap should beat or match the median random sample
  expect_gte(min_dist_wide, median(min_dists_random))
})

test_that("maximin with pre-computed dist works", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:29, nTip = 8)
  d <- TreeDist::ClusteringInfoDistance(trees)

  set.seed(6472)
  result_fn <- WideSample(trees, 5)
  set.seed(6472)
  result_dist <- WideSample(trees, 5, distance = d)

  expect_equal(result_fn, result_dist)
})

test_that("attributes are preserved", {
  trees <- as.phylo(0:9, nTip = 8)
  attr(trees, "score") <- 123.4
  attr(trees, "hits_to_best") <- 5L
  set.seed(2847)
  result <- WideSample(trees, 3)
  expect_equal(attr(result, "score"), 123.4)
  expect_equal(attr(result, "hits_to_best"), 5L)
})

test_that("input validation works", {
  trees <- as.phylo(0:4, nTip = 8)
  expect_error(WideSample(list(1, 2), 2), "multiPhylo")
  expect_error(WideSample(trees, -1), "non-negative")
  expect_error(WideSample(trees, c(1, 2)), "single")
  expect_error(WideSample(trees, NA), "non-negative")
})

test_that("dist dimension mismatch is caught", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:9, nTip = 8)
  bad_dist <- TreeDist::ClusteringInfoDistance(trees[1:5])
  expect_error(WideSample(trees, 3, distance = bad_dist), "entries")
})

test_that("bad distance argument is caught", {
  trees <- as.phylo(0:4, nTip = 8)
  expect_error(WideSample(trees, 2, distance = "nope"), "function or a dist")
})

test_that("maximin is reproducible with set.seed", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:49, nTip = 10)
  set.seed(3319)
  r1 <- WideSample(trees, 8)
  set.seed(3319)
  r2 <- WideSample(trees, 8)
  expect_equal(r1, r2)
})
