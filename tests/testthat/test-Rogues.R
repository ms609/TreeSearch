test_that("Rogues found with TipInstability()", {
  
  library("TreeTools", quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  instab <- TipInstability(trees)
  expect_equal('Rogue', names(which.max(instab)))
})
