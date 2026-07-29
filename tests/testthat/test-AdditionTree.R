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

  # `sequence` defaults to a random addition order, so seed for a reproducible
  # tree shape rather than for a lucky one: under T-364 these two assertions were
  # order-dependent and this seed was load-bearing, which is why the sweep in the
  # next test -- not this seed -- is now the guard.
  set.seed(1)
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

test_that("Addition tree obeys constraints for every addition order", {
  # T-364/T-370.  AdditionTree() seeds its search with a three-taxon tree built
  # from the first three taxa of the addition order, before any constraint is
  # consulted.  When that seed puts a constrained group's taxa on both sides of
  # its root, the group's LCA is the root itself; the constraint used to be
  # skipped from that point on -- and an LCA never moves back down, so it stayed
  # skipped for every remaining insertion and the tree came back violating the
  # constraint with no warning.  Enforcing through the split's complement fixes
  # it, since the same unrooted bipartition is displayed either way.
  #
  # Sweeps rather than single calls, because the failure was order-dependent: on
  # the code before the fix these three assertions report 35/400 seeds, 16/120
  # base triples and 42/120 base triples respectively, so a single unseeded call
  # passed ~91% of the time and a pinned seed proved nothing about any other
  # platform.
  dataset <- TreeTools::MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  taxa <- letters[1:6]
  efSplit <- as.Splits(c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE), taxa)
  efConstraint <- TreeTools::MatrixToPhyDat(
    c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1))

  # Random addition orders: seeds 23, 25, 29, 43 and 91 were the first that put
  # e and f in the base tree astride its root.
  broken <- which(vapply(seq_len(400), function(seed) {
    set.seed(seed)
    !(efSplit %in% as.Splits(AdditionTree(dataset, constraint = efConstraint),
                             taxa))
  }, logical(1)))
  expect_equal(broken, integer(0))

  # Every arrangement of the base triple, via `sequence`.  A group of three
  # covers the case where the triple falls entirely inside the constrained
  # group, which no rearrangement of those three taxa can rescue: whichever of
  # them the base tree puts opposite the root is itself inside the group, so the
  # group can only be enforced through its complement.
  defConstraint <- TreeTools::MatrixToPhyDat(
    c(a = 0, b = 0, c = 0, d = 1, e = 1, f = 1))
  defSplit <- as.Splits(c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), taxa)
  triples <- expand.grid(taxa, taxa, taxa, stringsAsFactors = FALSE)
  triples <- triples[apply(triples, 1, anyDuplicated) == 0L, ]
  BadOrders <- function(cons, split) {
    bad <- vapply(seq_len(nrow(triples)), function(i) {
      triple <- unlist(triples[i, ], use.names = FALSE)
      order <- c(triple, setdiff(taxa, triple))
      tree <- AdditionTree(dataset, constraint = cons, sequence = order)
      if (split %in% as.Splits(tree, taxa)) "" else paste(order, collapse = "")
    }, character(1))
    bad[nzchar(bad)]
  }
  expect_equal(BadOrders(efConstraint, efSplit), character(0))
  expect_equal(BadOrders(defConstraint, defSplit), character(0))

  # Honouring the constraint must not cost a taxon every legal insertion edge:
  # exhausting them falls back to an unchecked edge, which warns.  (This passed
  # before the fix too -- the old failure was silent -- so it guards against the
  # fix over-constraining, not against T-364 itself.)
  set.seed(23)
  expect_no_warning(AdditionTree(dataset, constraint = efConstraint))
  expect_no_warning(
    AdditionTree(dataset, constraint = defConstraint,
                 sequence = c("d", "e", "f", "a", "b", "c"))
  )
})

test_that("AdditionTree() rooting is an arbitrary construction artefact", {
  library("TreeTools", quietly = TRUE)
  # 6-taxon dataset from a Fitch phylogeny; sequence[1] is NOT reliably
  # the root -- documenting @return should not promise otherwise.
  dataset <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1,
      1, 0, 0, 1, 1, 0), ncol = 3,
    dimnames = list(letters[1:6], NULL)))

  RootTipLabel <- function (tr) {
    rootNode <- RootNode(tr)
    kids <- tr$edge[tr$edge[, 1] == rootNode, 2]
    tipKids <- kids[kids <= NTip(tr)]
    if (length(tipKids)) tr$tip.label[tipKids] else NA_character_
  }

  set.seed(1)
  tr <- AdditionTree(dataset, sequence = letters[1:6])
  # sequence[1] ("a") should not be assumed to be the root tip
  expect_false(identical(RootTipLabel(tr), "a"))
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

test_that("AdditionTree() verifies its own output against the constraint", {
  # The placement filter is not trusted to have been exhaustive: the finished
  # tree is checked against every constraint split and a warning raised if any
  # is missing.  `constraint_fallback` alone never sufficed -- it only fires
  # when the filter rejected *every* edge, which the T-364/T-370 leak never did,
  # so violating trees came back mutely -- and AdditionTree() never sets
  # `has_posthoc`, so unlike the search path it has no reshuffle to fall back
  # on.  Disabling complement enforcement makes 425 of 1334 randomised cases
  # violate; the check caught all 425 and warned on none of the other 909.
  #
  # The risk of adding a verifier is spurious warnings, so sweep for silence.
  dataset <- TreeTools::MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  efConstraint <- TreeTools::MatrixToPhyDat(
    c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1))

  warnings <- capture_warnings(
    for (seed in seq_len(150)) {
      set.seed(seed)
      AdditionTree(dataset, constraint = efConstraint)
    }
  )
  expect_equal(warnings, character(0))
})
