test_that("Maximum entropy calculated ok", {
  expect_equal(.MaxEntropy(5), 0)
  expect_equal(.MaxEntropy(rep(c("A", "B"), c(3, 5))), .MaxEntropy(c(5, 3)))
  expect_equal(.MaxEntropy(c(4, 4)), .Entropy(c(2, 2, 2, 2)))
  expect_equal(.MaxEntropy(c(5, 3)), .Entropy(c(3, 2, 2, 1)))
  expect_equal(.MaxEntropy(c(7, 3)), .Entropy(c(4, 3, 2, 1)))
  expect_equal(.MaxEntropy(c(2, 8)), .Entropy(c(4, 2, 4)))
})
