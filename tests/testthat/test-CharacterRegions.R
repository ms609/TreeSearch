test_that("CharacterRegions() behaves", {
  tokens <- upper.tri(matrix(0, 6, 6)) * 1
  rownames(tokens) <- TipLabels(6)
  mat <- MatrixToPhyDat(tokens)
  CharacterRegions(TreeTools::BalancedTree(6), mat)
})