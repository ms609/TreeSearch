test_that("SPR errors", {
  expect_error(asan_error(matrix(9, 1, 1)))
})