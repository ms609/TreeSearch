test_that("CharacterRegions() behaves", {
library("TreeTools") #TODO DELETE
  tokens <- cbind(upper.tri(matrix(0, 6, 6)) * 1, c(0, 1, 1, 1, 1, 0))
  rownames(tokens) <- TipLabels(6)
  mat <- MatrixToPhyDat(tokens)
  tr <- TreeTools::BalancedTree(6)
  tr$edge
  if (!interactive()) cr <- CharacterRegions(tr, mat)
  expect_equal(length(cr), ncol(tokens))
  expect_equal(cr, list(6, c(5, 1), c(4, 2), c(3, 3), c(3, 2, 1), c(5, 1), c(4, 1, 1)))
    
})
