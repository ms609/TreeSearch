# Tier 1: runs on CRAN
# Tests for WideSample() — Max-Min diversity (MMDP) tree subsampling

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

test_that("n = 1 returns a single tree deterministically", {
  trees <- as.phylo(0:9, nTip = 8)
  result <- WideSample(trees, 1)
  expect_equal(length(result), 1)
  expect_equal(result, WideSample(trees, 1))  # deterministic, no RNG
})

test_that("WideSample selects more diverse trees than random draws", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:99, nTip = 12)
  n <- 10
  # Minimum pairwise distance within a tree subset (the T_k objective).
  minGap <- function(treeSub) {
    m <- as.matrix(TreeDist::ClusteringInfoDistance(treeSub))
    min(m[lower.tri(m)])
  }

  wideGap <- minGap(WideSample(trees, n, quality = 1))

  set.seed(5028)
  randomGaps <- vapply(1:20, function(i)
    minGap(trees[sample.int(length(trees), n)]), double(1))
  # The dispersion-maximising subset should beat the median random draw.
  expect_gte(wideGap, median(randomGaps))
})

test_that("pre-computed dist (object and matrix) match the function path", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:29, nTip = 8)
  d <- TreeDist::ClusteringInfoDistance(trees)

  result_fn   <- WideSample(trees, 5, quality = 1)
  result_dist <- WideSample(trees, 5, dist = d, quality = 1)
  result_mat  <- WideSample(trees, 5, dist = as.matrix(d), quality = 1)

  expect_equal(result_fn, result_dist)
  expect_equal(result_fn, result_mat)
})

test_that("attributes are preserved", {
  trees <- as.phylo(0:9, nTip = 8)
  attr(trees, "score") <- 123.4
  attr(trees, "hits_to_best") <- 5L
  result <- WideSample(trees, 3)
  expect_equal(attr(result, "score"), 123.4)
  expect_equal(attr(result, "hits_to_best"), 5L)
})

test_that("WideSample is deterministic (no RNG)", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:49, nTip = 10)
  expect_equal(WideSample(trees, 8), WideSample(trees, 8))
})

# Input validation --------------------------------------------------------

test_that("input validation works", {
  trees <- as.phylo(0:4, nTip = 8)
  expect_error(WideSample(list(1, 2), 2), "multiPhylo")
  expect_error(WideSample(trees, -1), "non-negative")
  expect_error(WideSample(trees, c(1, 2)), "single")
  expect_error(WideSample(trees, NA), "non-negative")
})

test_that("quality is validated", {
  trees <- as.phylo(0:9, nTip = 8)
  expect_error(WideSample(trees, 3, quality = 4), "quality")
  expect_error(WideSample(trees, 3, quality = 0), "quality")
})

test_that("dist dimension mismatch is caught", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:9, nTip = 8)
  bad_dist <- TreeDist::ClusteringInfoDistance(trees[1:5])
  expect_error(WideSample(trees, 3, dist = bad_dist), "rows")
})

test_that("bad dist argument is caught", {
  trees <- as.phylo(0:4, nTip = 8)
  expect_error(WideSample(trees, 2, dist = "nope"),
               "must be a function")
})

# Solver tiers ------------------------------------------------------------

test_that("quality 1/2 return valid diverse subsets", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:39, nTip = 10)
  for (q in 1:2) {
    res <- WideSample(trees, 6, quality = q, time_budget_s = 1)
    expect_equal(length(res), 6)
    expect_true(all(names(res) %in% names(trees)))
  }
})

test_that("quality 3 matches or beats the heuristic optimum on small sets", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("highs")
  trees <- as.phylo(0:14, nTip = 8)
  n <- 4
  minGap <- function(treeSub) {
    m <- as.matrix(TreeDist::ClusteringInfoDistance(treeSub))
    min(m[lower.tri(m)])
  }
  exact     <- WideSample(trees, n, quality = 3)
  heuristic <- WideSample(trees, n, quality = 1)
  expect_gte(minGap(exact), minGap(heuristic))  # exact optimum is best-or-equal
})

test_that("forced quality 3 on a large set warns about cost", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("highs")
  skip_on_cran()
  trees <- as.phylo(0:44, nTip = 8)  # > 40 trips the warning
  expect_warning(
    WideSample(trees, 3, quality = 3, time_budget_s = 1),
    "may be very slow"
  )
})

# Tier-selection logic (no solving) ---------------------------------------

test_that(".SelectWideSampleTier keys on (matrix-available, N)", {
  sel <- TreeSearch:::.SelectWideSampleTier
  ceil <- 12000L
  # Auto: matrix present beyond the ceiling still reaches DropAdd, not Gonzalez.
  expect_equal(sel(NULL, TRUE,  20000L, ceil), 2L)
  # Auto: a function beyond the ceiling falls back to anchored Gonzalez.
  expect_equal(sel(NULL, FALSE, 20000L, ceil), 1L)
  # Auto: moderate N builds the matrix for DropAdd.
  expect_equal(sel(NULL, FALSE, 5000L,  ceil), 2L)
  # Forced quality 1 is always reachable.
  expect_equal(sel(1L,   FALSE, 20000L, ceil), 1L)
  # Forced quality 2/3 with a function past the ceiling is a hard error.
  expect_error(sel(2L, FALSE, 20000L, ceil), "build ceiling")
  expect_error(sel(3L, FALSE, 20000L, ceil), "build ceiling")
  # Forced quality 2/3 with a supplied matrix past the ceiling is fine.
  expect_equal(sel(2L, TRUE, 20000L, ceil), 2L)
  expect_equal(sel(3L, TRUE, 20000L, ceil), 3L)
})

test_that(".WideSampleColumnOracle probes the (tree, trees) contract", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:9, nTip = 8)
  oracle <- TreeSearch:::.WideSampleColumnOracle(
    TreeDist::ClusteringInfoDistance, trees, length(trees)
  )
  col1 <- oracle(1L)
  expect_length(col1, length(trees))
  expect_true(is.numeric(col1))

  # A function that returns the wrong length is rejected up front.
  badFn <- function(a, b) 1:3
  expect_error(
    TreeSearch:::.WideSampleColumnOracle(badFn, trees, length(trees)),
    "length"
  )
})
