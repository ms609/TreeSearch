context("combinatorics.R")

test_that("Factorials are calculated correctly", {
  expect_equal(c(1L, 1L, 1L, 2L, 3L, 2L * 4L, 3L * 5L,
                 2L * 4L * 6L, 3L * 5L * 7L), DoubleFactorial(-1:7))
  expect_equal(doubleFactorials[1:20], DoubleFactorial(1:20))
  expect_equal(LogDoubleFactorial(-1:10), log(DoubleFactorial(-1:10)))
  expect_equal(logDoubleFactorials[1:20], logDoubleFactorials[1:20])
})