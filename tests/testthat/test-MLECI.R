test_that("MLCI works", {
  tree <- TreeTools::BalancedTree(8)
  chars <- cbind(c(0, 0, 0, 0, 1, 1, 1, 1),
                 c(0, 0, 0, 0, 1, 0, 0, 0),
                 c(0, 0, 0, 0, 0, 1, 1, 1),
                 c(0, 1, 0, 1, 0, 1, 0, 1))
  
  expect_no_error(val <- MLCI(tree, chars, precision = 1e-2))
  val
  expect_equal(val[1, "mlci"], 1)
  expect_equal(val[2, "mlci"], 0)
  expect_gt(val[1, "mlci"], val[3, "mlci"])
  expect_lt(val[4, "mlci"], 0)
})
