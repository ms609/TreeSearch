test_that("MLCI works", {
  tree <- TreeTools::BalancedTree(8)
  chars <- cbind(c(0, 0, 0, 0, 1, 1, 1, 1),
                c(0, 0, 0, 0, 1, 0, 0, 0),
                c(0, 0, 0, 0, 0, 1, 1, 1),
                c(0, 1, 0, 1, 0, 1, 0, 1))
  expect_no_error(MLCI(tree, chars))
})
