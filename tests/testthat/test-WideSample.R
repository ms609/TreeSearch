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

test_that("n = 1 returns the medoid deterministically", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:9, nTip = 8)
  names(trees) <- paste0("t", seq_along(trees))
  d <- as.matrix(TreeDist::ClusteringInfoDistance(trees))
  result <- WideSample(trees, 1, dist = d)
  expect_equal(length(result), 1)
  # The medoid minimizes summed distance to the rest (most central tree).
  expect_equal(names(result), names(trees)[which.min(rowSums(d))])
  expect_equal(result, WideSample(trees, 1, dist = d))  # deterministic, no RNG
})

test_that("n = 1 medoid falls back to a matrix-free seed past the ceiling", {
  skip_if_not_installed("TreeDist")
  old <- options(TreeSearch.WideSample.buildCeiling = 5L)
  on.exit(options(old))
  trees <- as.phylo(0:9, nTip = 8)  # 10 > the (lowered) build ceiling of 5
  names(trees) <- paste0("t", seq_along(trees))
  result <- WideSample(trees, 1)    # function path, no matrix built
  expect_equal(length(result), 1)
  expect_true(all(names(result) %in% names(trees)))
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

  wideGap <- minGap(WideSample(trees, n, effort = 1))

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

  result_fn   <- WideSample(trees, 5, effort = 1)
  result_dist <- WideSample(trees, 5, dist = d, effort = 1)
  result_mat  <- WideSample(trees, 5, dist = as.matrix(d), effort = 1)

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

test_that("WideSample is deterministic on the RNG-free tiers", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:49, nTip = 10)
  # FarFirst (peripheral seed) and DropAdd (tabu search) draw on no RNG, so a
  # forced effort 1 or 2 is bit-reproducible. (Auto, Grasp and exact may vary
  # their selection with the session RNG -- see the Grasp test below.)
  expect_equal(WideSample(trees, 8, effort = 1),
               WideSample(trees, 8, effort = 1))
  expect_equal(WideSample(trees, 8, effort = 2),
               WideSample(trees, 8, effort = 2))
})

test_that("Grasp tier (effort 3) is reproducible under set.seed", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:39, nTip = 10)
  set.seed(1)
  a <- WideSample(trees, 6, effort = 3)
  set.seed(1)
  b <- WideSample(trees, 6, effort = 3)
  expect_equal(a, b)
})

# Input validation --------------------------------------------------------

test_that("input validation works", {
  trees <- as.phylo(0:4, nTip = 8)
  expect_error(WideSample(list(1, 2), 2), "multiPhylo")
  expect_error(WideSample(trees, -1), "non-negative")
  expect_error(WideSample(trees, c(1, 2)), "single")
  expect_error(WideSample(trees, NA), "non-negative")
})

test_that("effort is validated", {
  trees <- as.phylo(0:9, nTip = 8)
  expect_error(WideSample(trees, 3, effort = 5), "effort")
  expect_error(WideSample(trees, 3, effort = 0), "effort")
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

test_that("effort 1/2/3 return valid diverse subsets", {
  skip_if_not_installed("TreeDist")
  trees <- as.phylo(0:39, nTip = 10)
  set.seed(1)  # effort 3 (Grasp) draws on the session RNG
  for (eff in 1:3) {
    res <- WideSample(trees, 6, effort = eff, maxSeconds = 1)
    expect_equal(length(res), 6)
    expect_true(all(names(res) %in% names(trees)))
  }
})

test_that("effort 4 matches or beats the heuristic optimum on small sets", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("highs")
  trees <- as.phylo(0:14, nTip = 8)
  n <- 4
  minGap <- function(treeSub) {
    m <- as.matrix(TreeDist::ClusteringInfoDistance(treeSub))
    min(m[lower.tri(m)])
  }
  exact     <- WideSample(trees, n, effort = 4)
  heuristic <- WideSample(trees, n, effort = 1)
  expect_gte(minGap(exact), minGap(heuristic))  # exact optimum is best-or-equal
})

test_that("forced effort 4 above the exact ceiling warns about cost", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("highs")
  skip_on_cran()
  trees <- as.phylo(0:209, nTip = 8)  # > 200 trips the warning
  expect_warning(
    WideSample(trees, 3, effort = 4, maxSeconds = 1),
    "may be very slow"
  )
})

# Tier-selection logic (no solving) ---------------------------------------

test_that(".SelectWideSampleTier keys on (matrix-available, N)", {
  sel <- TreeSearch:::.SelectWideSampleTier
  ceil <- 12000L
  exCeil <- 200L
  # Auto: matrix present beyond the build ceiling still reaches DropAdd, not
  # the matrix-free FarFirst.
  expect_equal(sel(NULL, TRUE,  20000L, ceil), 2L)
  # Auto: a function beyond the build ceiling falls back to FarFirst.
  expect_equal(sel(NULL, FALSE, 20000L, ceil), 1L)
  # Auto: moderate N builds the matrix for DropAdd.
  expect_equal(sel(NULL, FALSE, 5000L,  ceil), 2L)
  # Auto: small N reaches the exact tier when highs is available, else DropAdd.
  expect_equal(sel(NULL, FALSE, 50L, ceil, exCeil, highsAvailable = TRUE), 4L)
  expect_equal(sel(NULL, FALSE, 50L, ceil, exCeil, highsAvailable = FALSE), 2L)
  # Auto: just above the exact ceiling never auto-selects exact (Grasp is
  # opt-in, so DropAdd is the auto heuristic here).
  expect_equal(sel(NULL, TRUE, 201L, ceil, exCeil, highsAvailable = TRUE), 2L)
  # Forced effort 1 is always reachable.
  expect_equal(sel(1L,   FALSE, 20000L, ceil), 1L)
  # Forced effort 2/3/4 with a function past the build ceiling is a hard error.
  expect_error(sel(2L, FALSE, 20000L, ceil), "build ceiling")
  expect_error(sel(3L, FALSE, 20000L, ceil), "build ceiling")
  expect_error(sel(4L, FALSE, 20000L, ceil), "build ceiling")
  # Forced effort 2/3/4 with a supplied matrix past the ceiling is fine.
  expect_equal(sel(2L, TRUE, 20000L, ceil), 2L)
  expect_equal(sel(3L, TRUE, 20000L, ceil), 3L)
  expect_equal(sel(4L, TRUE, 20000L, ceil), 4L)
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
