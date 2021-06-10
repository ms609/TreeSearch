test_that("Rogues found with TipInstability()", {
  
  library("TreeTools", quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  instab <- TipInstability(trees)
  expect_equal('Rogue', names(which.max(instab)))
  mi <- TipInformation(trees)
  expect_equal('Rogue', names(mi[mi > 0]))
  
})

test_that("Wilkinson & Crotti's examples are satisfied", {
  fig2 <- list(AddTip(BalancedTree(6:1), '3', 'X'),
               AddTip(BalancedTree(6:1), '4', 'X'))
  
})
