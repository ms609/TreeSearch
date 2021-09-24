test_that("SPR errors", {
  library("TreeTools")
  expect_error(asan_error(as.phylo(1, 3)$edge)), integer(0)
  expect_error(asan_error(Postorder(as.phylo(1, 6))$edge, integer(0)))
  expect_error(asan_error(SortTree(as.phylo(1, 6))$edge, integer(0)))
})