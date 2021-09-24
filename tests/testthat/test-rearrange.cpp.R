test_that("SPR errors", {
  asan_error(matrix(9, 1, 1))
})