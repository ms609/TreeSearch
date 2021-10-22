test_that("ClusterStrings() works", {
  expect_equal(ClusterStrings(rep(letters[1:6], 1:6)),
               structure(rep(1:6, 1:6), 'med' = letters[1:6]))
  expect_equal(range(ClusterStrings(rep(letters[1:6], 1:6), 2)),
               1:2)
  expect_equal(ClusterStrings(paste0(c('aaaa', 'bbb', 'cccccc'), 1:20)),
               structure(rep_len(1:3, 20),
                         silhouette = 0.7955785, # copied, not calculated
                         med = paste0(c('aaaa', 'bbb', 'cccccc'), 1:3)))
})