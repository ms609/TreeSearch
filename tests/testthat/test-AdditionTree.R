test_that("Addition tree is more parsimonious", {
  data('Lobo', package = 'TreeTools')
  L10 <- Lobo.phy[1:10]
  seq10 <- names(L10)
  Score10 <- function (tr) TreeLength(tr, Lobo.phy, concavity = 10)
  
  eq <- AdditionTree(Lobo.phy)
  kx <- AdditionTree(L10, seq10, concavity = 10)
  #pr <- AdditionTree(Lobo.phy, concavity = 'pr')
  nj <- NJTree(Lobo.phy)
  
  expect_lt(TreeLength(eq, Lobo.phy), TreeLength(nj, Lobo.phy))
  expect_lt(TreeLength(kx, L10, concavity = 10),
            TreeLength(TreeTools::KeepTip(nj, 1:10), L10, concavity = 10))
  
  # expect_lt(TreeLength(pr, Lobo.phy, concavity = 'pr'),
  #           TreeLength(nj, Lobo.phy, concavity = 'pr'))
  
})
