test_that("Addition tree produces valid trees", {
  data("Lobo", package = "TreeTools")
  L10 <- Lobo.phy[1:10]
  seq10 <- names(L10)

  set.seed(1)
  eq <- AdditionTree(Lobo.phy)
  expect_equal(TreeTools::NTip(eq), length(Lobo.phy))
  expect_true(TreeLength(eq, Lobo.phy) > 0)

  kx <- AdditionTree(L10, sequence = seq10, concavity = 10)
  expect_equal(TreeTools::NTip(kx), 10L)

  # PrepareDataProfile() emits a cli message about inapplicable tokens for
  # profile parsimony; suppress so it doesn't leak into testthat output.
  pr <- suppressMessages(
    AdditionTree(L10, sequence = 1:10, concavity = "profile")
  )
  expect_equal(TreeTools::NTip(pr), 10L)
})

test_that(".ConstraintConstrains() succeeds", {
  expect_false(TreeSearch:::.ConstraintConstrains(NULL))

  # Single level
  expect_false(TreeSearch:::.ConstraintConstrains(
    structure(list(A = 1L, B = 2L, C = 2L, D = 2L), weight = 1L, nr = 1L,
              nc = 1L, index = 1L, levels = 0, allLevels = c("0", "?"),
              type = "USER", contrast =
                structure(c(1, 1), dim = 2:1, dimnames = list(NULL, 0)),
              class = "phyDat")
  ))

  expect_false(TreeSearch:::.ConstraintConstrains(
    structure(list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L), weight = 1L, nr = 1L,
              nc = 2L, index = 1L, levels = 0:1,
              allLevels = c("0", "1", "?"), type = "USER",
              contrast = structure(c(1, 0, 1, 0, 1, 1), dim = 3:2,
                                   dimnames = list(NULL, 0:1)),
              class = "phyDat")
  ))
  expect_true(TreeSearch:::.ConstraintConstrains(structure(
    list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L, F = 2L), weight = 1L, nr = 1L,
    nc = 2L, index = 1L, levels = 0:1, allLevels = c("0", "1", "?"),
    type = "USER", contrast = structure(c(1, 0, 1, 0, 1, 1), dim = 3:2,
                                        dimnames = list(NULL, 0:1)),
    class = "phyDat")
  ))
  expect_false(TreeSearch:::.ConstraintConstrains(structure(
    list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L, F = 2L), weight = 1L, nr = 1L,
    nc = 2L, index = 1L,
    levels = 0:2, allLevels = c("0", "1", "2", "?"), type = "USER",
    contrast = structure(c(1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1),
                         dim = c(4, 3), dimnames = list(NULL, 0:2)),
    class = "phyDat")
  ))
  expect_true(TreeSearch:::.ConstraintConstrains(structure(
    list(A = 1L, B = 2L, C = 1L, D = 1L, E = 3L, F = 2L), weight = 1L, nr = 1L,
    nc = 2L, index = 1L, levels = 0:2, allLevels = c("0", "1", "2", "?"),
    type = "USER", contrast = structure(c(1, 0, 1, 1, 0, 1, 0, 1, 1),
                                        dim = c(3, 3), dimnames = list(NULL, 0:2)),
    class = "phyDat")
  ))
})

test_that("Addition tree obeys constraints", {
  dataset <- TreeTools::MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  constraint <- c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1)
  expected_split <- as.Splits(c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
                               letters[1:6])

  # as phyDat
  expect_true(expected_split %in%
              as.Splits(AdditionTree(dataset,
                constraint = TreeTools::MatrixToPhyDat(constraint)),
                letters[1:6]))
  # as matrix
  expect_true(expected_split %in%
              as.Splits(AdditionTree(dataset, constraint = cbind(constraint)),
                letters[1:6]))

  # Trivial constraints should not affect tree
  set.seed(0)
  unconstrained <- AdditionTree(dataset)

  set.seed(0)
  expect_equal(AdditionTree(dataset, constraint = NULL), unconstrained)

  # Partial constraint with subset of taxa
  cdef <- letters[3:6]
  set.seed(0)
  subtree <- TreeTools::KeepTip(
    AdditionTree(dataset, constraint = constraint[3:6], seq = letters[1:6]),
    cdef)
  expect_equal_tree(ape::read.tree(text = "(c, d, (e, f));"),
               TreeTools::UnrootTree(subtree))
})

test_that("AdditionTree() handles edge cases", {
  library("TreeTools", quietly = TRUE)
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  expect_equal(PectinateTree(letters[1:3]), AdditionTree(dataset[1:3]))
  expect_equal(5, NTip(AdditionTree(dataset[-4])))
  # 4-tip profile tree — suppress PrepareDataProfile() cli message.
  expect_equal(4L, NTip(suppressMessages(
    AdditionTree(dataset[1:4], conc = "profile")
  )))
})

test_that("AdditionTree() rejects duplicated `sequence` taxa", {
  library("TreeTools", quietly = TRUE)
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  taxa <- names(dataset)

  # A duplicated taxon name in a *character* sequence used to slip past
  # validation and poison the C++ kernel's addition order: the repeated tip
  # was inserted twice and a different tip never added, so AdditionTree()
  # silently returned a phylo containing one taxon twice and dropping another
  # (which still passed checkValidPhylo / is.binary).  The numeric path always
  # rejected duplicates; the character path must too.
  expect_error(AdditionTree(dataset, sequence = c(taxa[1], taxa[1], taxa[2])),
               "more than once")
  expect_error(
    AdditionTree(dataset, sequence = c(taxa[1], taxa[2:5], taxa[1])),
    "more than once")
  # numeric duplicates remain rejected (regression guard for both paths)
  expect_error(AdditionTree(dataset, sequence = c(1L, 1L, 2L)),
               "distinct whole-number")

  # Valid distinct sequences (full + partial) are unaffected.
  expect_equal(NTip(AdditionTree(dataset, sequence = taxa)), 6L)
  expect_equal(NTip(AdditionTree(dataset, sequence = taxa[c(3, 1)])), 6L)
})
