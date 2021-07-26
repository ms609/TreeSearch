test_that("Addition tree is more parsimonious", {
  data('Lobo', package = 'TreeTools')
  L10 <- Lobo.phy[1:10]
  seq10 <- names(L10)
  Score10 <- function (tr) TreeLength(tr, Lobo.phy, concavity = 10)
  
  set.seed(1) # ensure consistent addition sequence
  eq <- AdditionTree(Lobo.phy)
  kx <- AdditionTree(L10, sequence = seq10, concavity = 10)
  #pr <- AdditionTree(Lobo.phy, concavity = 'pr')
  nj <- TreeTools::NJTree(Lobo.phy)
  
  expect_lt(TreeLength(eq, Lobo.phy), TreeLength(nj, Lobo.phy))
  expect_lt(TreeLength(kx, L10, concavity = 10),
            TreeLength(TreeTools::KeepTip(nj, 1:10), L10, conc = 10))
  
  # expect_lt(TreeLength(pr, Lobo.phy, concavity = 'pr'),
  #           TreeLength(nj, Lobo.phy, concavity = 'pr'))
  
})

test_that("Addition tree obeys constraints", {
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  constraint <- MatrixToPhyDat(c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1))
  expect_true(as.Splits(c(F, F, F, F, T, T), letters[1:6]) %in% 
              as.Splits(AdditionTree(dataset, constraint = constraint),
                        letters[1:6]))
  
  cdef <- letters[3:6]
  subtree <- TreeTools::KeepTip(
    AdditionTree(dataset, constraint = constraint[3:6], seq = letters[1:6]), 
    cdef)
  expect_equal(read.tree(text = '(c, d, (e, f));'),
               TreeTools::UnrootTree(subtree))
})
