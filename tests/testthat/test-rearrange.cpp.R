test_that("SPR errors", {
  library("TreeTools")
  expect_error(asan_error(as.phylo(1, 3)$edge))
})