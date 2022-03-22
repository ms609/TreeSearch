test_that("ClusterStrings() works", {
  skip("Locate memory leak, #TODO remove")
  x <- rep(letters[1:6], 1:6)
  expect_equal(ClusterStrings(x),
               structure(rep(1:6, 1:6), 'med' = letters[1:6]))
  expect_error(ClusterStrings(x, 1), "`maxCluster` must be at least two.")
  expect_equal(range(ClusterStrings(x, 2)), 1:2)
  expect_equal(ClusterStrings(paste0(c('aaaa', 'bbb', 'cccccc'), 1:20)),
               structure(rep_len(1:3, 20),
                         silhouette = 0.7955785, # copied, not calculated
                         med = paste0(c('aaaa', 'bbb', 'cccccc'), 1:3)))
})