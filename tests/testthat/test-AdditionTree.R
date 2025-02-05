test_that("Addition tree is more parsimonious", {
  data("Lobo", package = "TreeTools")
  L10 <- Lobo.phy[1:10]
  seq10 <- names(L10)
  Score <- function (tr, k) TreeLength(tr, Lobo.phy, concavity = k)
  
  set.seed(1) # ensure consistent addition sequence
  eq <- AdditionTree(Lobo.phy)
  kx <- AdditionTree(L10, sequence = seq10, concavity = 10)
  pr <- AdditionTree(L10, sequence = 1:10, concavity = "pr")
  
  skip_if_not_installed("phangorn")
  # Previously used TreeTools::NJTree but since rewriting it's more parsimonious
  # than ape/phangorn.
  nj <- RootTree(ape::nj(phangorn::dist.hamming(Lobo.phy)), 1)
  nj10 <- TreeTools::KeepTip(nj, 1:10)
  
  expect_lt(TreeLength(eq, Lobo.phy), TreeLength(nj, Lobo.phy))
  expect_lt(Score(kx, 10), Score(nj10, 10))
  expect_lt(Score(pr, "pr"), Score(nj10, "pr"))
})

test_that(".ConstraintConstrains() succeeds", {
  expect_false(.ConstraintConstrains(NULL))
  
  # Single level
  expect_false(.ConstraintConstrains(
    structure(list(A = 1L, B = 2L, C = 2L, D = 2L), weight = 1L, nr = 1L,
              nc = 1L, index = 1L, levels = 0, allLevels = c("0", "?"),
              type = "USER", contrast = 
                structure(c(1, 1), dim = 2:1, dimnames = list(NULL, 0)),
              class = "phyDat")
  ))
  
  expect_false(.ConstraintConstrains(
    structure(list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L), weight = 1L, nr = 1L,
              nc = 2L, index = 1L, levels = 0:1,
              allLevels = c("0", "1", "?"), type = "USER",
              contrast = structure(c(1, 0, 1, 0, 1, 1), dim = 3:2,
                                   dimnames = list(NULL, 0:1)),
              class = "phyDat")
  ))
  expect_true(.ConstraintConstrains(structure(
    list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L, F = 2L), weight = 1L, nr = 1L,
    nc = 2L, index = 1L, levels = 0:1, allLevels = c("0", "1", "?"),
    type = "USER", contrast = structure(c(1, 0, 1, 0, 1, 1), dim = 3:2,
                                        dimnames = list(NULL, 0:1)),
    class = "phyDat")
  ))
  expect_false(.ConstraintConstrains(structure(
    list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L, F = 2L), weight = 1L, nr = 1L,
    nc = 2L, index = 1L,
    levels = 0:2, allLevels = c("0", "1", "2", "?"), type = "USER",
    contrast = structure(c(1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1),
                         dim = c(4, 3), dimnames = list(NULL, 0:2)),
    class = "phyDat")
  ))
  expect_true(.ConstraintConstrains(structure(
    list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L, F = 2L), weight = 1L, nr = 1L,
    nc = 2L, index = 1L, levels = 0:2, allLevels = c("0", "1", "2", "?"),
    type = "USER", contrast = structure(c(1, 0, 1, 1, 0, 1, 0, 1, 1),
                                        dim = c(3, 3), dimnames = list(NULL, 0:2)),
    class = "phyDat")
  ))
})

test_that("Addition tree obeys constraints", {
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  constraint <- c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1)
  # as phyDat
  expect_true(as.Splits(c(F, F, F, F, T, T), letters[1:6]) %in%
              as.Splits(AdditionTree(dataset, constraint = 
                                       MatrixToPhyDat(constraint)),
                        letters[1:6]))
  # as non-phyDat
  expect_true(as.Splits(c(F, F, F, F, T, T), letters[1:6]) %in%
              as.Splits(AdditionTree(dataset, constraint = cbind(constraint)),
                        letters[1:6]))
  
  constraintTree <- TreeTools::BalancedTree(constraint)
  
  set.seed(0)
  unconstrained <- AdditionTree(dataset)
  
  CheckUnconstrained <- function(constraint) {
    set.seed(0)
    expect_equal(AdditionTree(dataset, constraint = constraint), unconstrained)
  }
  
  CheckUnconstrained(KeepTip(constraintTree, c("a", "b")))
  CheckUnconstrained(c(a = 0))
  CheckUnconstrained(KeepTip(constraintTree, "a"))
  CheckUnconstrained(c())
  CheckUnconstrained(KeepTip(constraintTree, character(0)))
  CheckUnconstrained(NULL)
  
  cdef <- letters[3:6]
  subtree <- TreeTools::KeepTip(
    AdditionTree(dataset, constraint = constraint[3:6], seq = letters[1:6]), 
    cdef)
  expect_equal(ape::read.tree(text = "(c, d, (e, f));"),
               TreeTools::UnrootTree(subtree))
})

test_that("AdditionTree() handles edge cases", {
  library("TreeTools", quietly = TRUE)
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  expect_equal(PectinateTree(letters[1:3]), AdditionTree(dataset[1:3]))
  expect_equal(UnrootTree(PectinateTree(c("a", "d", "b", "c"))), 
               UnrootTree(AdditionTree(dataset[1:4], conc = "pr")))
  # All trees have equal score
  expect_equal(5, NTip(AdditionTree(dataset[-4])))
})