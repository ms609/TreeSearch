test_that("Pedigree tree is more parsimonious", {
  data("Lobo", package = "TreeTools")
  L10 <- Lobo.phy[1:10]
  seq10 <- names(L10)
  Score <- function (tr, k) TreeLength(tr, Lobo.phy, concavity = k)
  
  set.seed(1) # ensure consistent sequence
  eq <- PedigreeTree(Lobo.phy)
  kx <- PedigreeTree(L10, sequence = seq10, concavity = 10)
  pr <- PedigreeTree(L10, sequence = 1:10, concavity = "pr")
  
  skip_if_not_installed("phangorn")
  # Previously used TreeTools::NJTree but since rewriting it's more parsimonious
  # than ape/phangorn.
  nj <- RootTree(ape::nj(phangorn::dist.hamming(Lobo.phy)), 1)
  nj10 <- TreeTools::KeepTip(nj, 1:10)
  
  expect_lt(TreeLength(eq, Lobo.phy), TreeLength(nj, Lobo.phy))
  expect_lt(Score(kx, 10), Score(nj10, 10))
  expect_lt(Score(pr, "pr"), Score(nj10, "pr"))
})

test_that("Pedigree tree obeys constraints", {
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  constraint <- c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1)
  # as phyDat
  expect_true(as.Splits(c(F, F, F, F, T, T), letters[1:6]) %in%
                as.Splits(PedigreeTree(dataset, constraint = 
                                         MatrixToPhyDat(constraint)),
                          letters[1:6]))
  # as non-phyDat
  expect_true(as.Splits(c(F, F, F, F, T, T), letters[1:6]) %in%
                as.Splits(PedigreeTree(dataset, constraint = cbind(constraint)),
                          letters[1:6]))
  
  cdef <- letters[3:6]
  subtree <- TreeTools::KeepTip(
    PedigreeTree(dataset, constraint = constraint[3:6], seq = letters[1:6]), 
    cdef)
  expect_equal(ape::read.tree(text = "(c, d, (e, f));"),
               TreeTools::UnrootTree(subtree))
})

test_that("PedigreeTree() handles edge cases", {
  library("TreeTools", quietly = TRUE)
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  expect_equal(PectinateTree(letters[1:3]), PedigreeTree(dataset[1:3]))
  expect_equal(UnrootTree(PectinateTree(c("a", "d", "b", "c"))), 
               UnrootTree(PedigreeTree(dataset[1:4], conc = "pr")))
  # All trees have equal score
  expect_equal(5, NTip(PedigreeTree(dataset[-4])))
})