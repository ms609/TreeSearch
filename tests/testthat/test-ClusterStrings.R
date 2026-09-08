skip_if_not_installed("cluster")
skip_if_not_installed("protoclust")

test_that("ClusterStrings() works", {
  x <- rep(letters[1:6], 1:6)
  expect_equal(ClusterStrings(x),
               structure(rep(1:6, 1:6), "med" = letters[1:6],
                         silhouette = NA_real_))
  expect_error(ClusterStrings(x, 1), "`maxCluster` must be at least two.")
  # Silhouette now computed on the true dissimilarity matrix (#97); the old
  # pam(dists, k) call clustered on Euclidean distance between rows of
  # `dists` instead, inflating this above the 0.5 "structure" threshold.
  expect_equal(range(ClusterStrings(x, 2)), c(1L, 1L))
  expect_equal(ClusterStrings(paste0(c("aaaa", "bbb", "cccccc"), 1:20)),
               structure(rep_len(1:3, 20),
                         # was 0.7955785 pre-fix
                         silhouette = 0.727540221,
                         med = paste0(c("aaaa", "bbb", "cccccc"), 1:3)),
               tolerance = 1e-6)
})

test_that("ClusterStrings() handles a singleton cluster (#96)", {
  # Pre-fix: colSums(dists[these, these]) drops to a scalar when the winning
  # clustering contains a singleton, erroring "'x' must be an array of at
  # least two dimensions".
  x <- c(paste0("aaaa", 1:5), paste0("bbbbbbbb", 1:5), paste0("cccccccccc", 1:5),
         "This is a totally different weird string zzzzzzzzzzzzzzzzzzzzzz")
  res <- expect_silent(ClusterStrings(x))
  expect_equal(length(res), length(x))
  expect_true(any(tabulate(res) == 1))
})

test_that("ClusterStrings() 'no structure' branch returns a per-element vector (#113)", {
  # Pre-fix: switch(..., 1) collapsed to a bare scalar `1` instead of
  # rep(1L, length(x)), violating the documented return contract.
  set.seed(42)
  x <- vapply(1:15, function(i) {
    paste(sample(letters, 6), collapse = "")
  }, character(1))
  res <- ClusterStrings(x)
  expect_equal(length(res), length(x))
  expect_true(all(res == 1L))
  expect_false(is.na(attr(res, "silhouette")))
})
