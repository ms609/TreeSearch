# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

## agent-issues/TreeSearch#54: the constraint the search enforces must
## be the one `?MaximizeParsimony`'s `constraint` argument documents.
##
## The documented contract is the phyDat reading: a returned tree is compliant
## when some edge separates the taxa coded `1` from those coded `0`, with
## `?`-coded taxa free to fall on either side.  The locked-node machinery used
## to enforce a strictly stronger one — some node's tip set had to EQUAL the 1
## group (or its exact complement), free taxa excluded.
##
## Strict implies loose, so no wrong answer was ever returned.  What broke was
## movement: a tree that satisfies the documented contract without making
## either group an exact clade mapped to no node at all, which
## regraft_violates_constraint() reads as "the tree already violates" and
## answers by rejecting EVERY rearrangement.  The replicate froze on its start.
##
## So these tests assert two things together, and neither alone is the fix:
## that the search MOVES from such a start, and that what it returns is still
## compliant.

library("TreeTools")

# 8 taxa; three characters agree on {a,b,e,f} | {c,d,g,h} and one cuts across
# it, so the start tree below is a local optimum only for a search that cannot
# move.
freeTaxaData <- function() {
  taxa <- letters[1:8]
  phangorn::phyDat(
    matrix(c("0", "0", "1", "1", "0", "0", "1", "1",
             "0", "0", "1", "1", "0", "0", "1", "1",
             "1", "1", "0", "0", "1", "1", "0", "0",
             "0", "1", "0", "1", "0", "1", "0", "1",
             "0", "1", "0", "1", "0", "1", "0", "1"),
           nrow = 8, dimnames = list(taxa, NULL)),
    type = "USER", levels = c("0", "1")
  )
}

# c(a = 1, b = 1, c = 0, d = 0, e:h = "?")
freeTaxaConstraint <- function() {
  phangorn::phyDat(
    matrix(c("1", "1", "0", "0", "?", "?", "?", "?"),
           nrow = 8, dimnames = list(letters[1:8], NULL)),
    type = "USER", levels = c("0", "1")
  )
}

# `{a,e,b,f}` | `{c,d,g,h}` separates {a,b} from {c,d}, so this satisfies the
# documented constraint — but {a,b} is not a clade, and neither is {c,d}.
freeTaxaStart <- function() {
  ape::read.tree(text = "(((a,e),(b,f)),(c,(d,(g,h))));")
}

# Does `tree` display a split with every tip of `inGroup` on one side and every
# tip of `outGroup` on the other?  Tips in neither group are ignored — this is
# the documented contract, spelled out independently of the engine.
SeparatesGroups <- function(tree, labels, inGroup, outGroup) {
  splits <- as.logical(as.Splits(tree, tipLabels = labels))
  if (!is.matrix(splits)) splits <- matrix(splits, nrow = 1)
  isIn <- labels %in% inGroup
  isOut <- labels %in% outGroup
  any(apply(splits, 1, function(row) {
    (all(row[isIn]) && !any(row[isOut])) ||
      (!any(row[isIn]) && all(row[isOut]))
  }))
}

# The engine hands back bare edge matrices for rooted binary trees.
EdgeSeparatesGroups <- function(edge, labels, inGroup, outGroup) {
  tree <- structure(
    list(edge = edge, Nnode = max(edge) - length(labels), tip.label = labels),
    class = "phylo")
  SeparatesGroups(tree, labels, inGroup, outGroup)
}

test_that("free `?` taxa do not freeze a compliant start tree", {
  dataset <- freeTaxaData()
  labels <- names(dataset)
  ds <- make_ts_data(dataset)
  start <- Preorder(RenumberTips(freeTaxaStart(), labels))
  startScore <- TreeLength(start, dataset)

  # Premise: the start satisfies the documented constraint but neither group is
  # a clade, so the pre-#54 exact-clade test could not map it.
  expect_true(SeparatesGroups(start, labels, c("a", "b"), c("c", "d")))
  expect_false(SeparatesGroups(start, labels, c("a", "b"),
                               setdiff(labels, c("a", "b"))))

  # Take the split matrix from .PrepareConstraint rather than writing it out,
  # so this exercises the same R -> C++ contract the user's `constraint =`
  # phyDat travels along: pre-#54 it coded the free taxa 0 (making {a,b} an
  # exact clade), now it codes them NA.
  free <- TreeSearch:::.PrepareConstraint(
    freeTaxaConstraint(), dataset)[["consSplitMatrix"]]
  expect_equal(as.vector(free), c(1L, 1L, 0L, 0L, NA, NA, NA, NA))

  set.seed(386)
  result <- tbrOnlyRun(ds, start[["edge"]], free)

  # The defect: every regraft was rejected and the start came back unimproved.
  expect_lt(result$best_score, startScore)

  # ... and what comes back still honours the documented constraint.
  expect_true(all(vapply(result$trees, EdgeSeparatesGroups, logical(1),
                         labels = labels, inGroup = c("a", "b"),
                         outGroup = c("c", "d"))))
})

test_that("a 0/1 constraint matrix still enforces the exact clade", {
  # Guard against over-loosening: with no free tips the two groups are
  # complements, and every check must collapse back to the exact-clade test.
  dataset <- freeTaxaData()
  labels <- names(dataset)
  ds <- make_ts_data(dataset)
  start <- Preorder(RenumberTips(
    ape::read.tree(text = "(((a,b),(e,f)),(c,(d,(g,h))));"), labels))

  strict <- matrix(c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L), nrow = 1)
  set.seed(386)
  result <- tbrOnlyRun(ds, start[["edge"]], strict)

  # Guard the guard: a search that froze would satisfy the compliance test
  # below for the wrong reason, since the start already has {a,b} as a clade.
  expect_lt(result$best_score, TreeLength(start, dataset))
  expect_true(all(vapply(result$trees, EdgeSeparatesGroups, logical(1),
                         labels = labels, inGroup = c("a", "b"),
                         outGroup = setdiff(labels, c("a", "b")))))
})

test_that("MaximizeParsimony honours and searches under a `?` constraint", {
  # End to end.  Unlike the TBR-only test above this cannot isolate the freeze
  # — ratchet, drift and nni-perturb all get a replicate moving again by other
  # means — so what it adds is that the whole pipeline still returns compliant
  # trees once every entry point reads the constraint the same way.
  dataset <- freeTaxaData()
  labels <- names(dataset)
  start <- freeTaxaStart()

  set.seed(386)
  result <- suppressWarnings(
    MaximizeParsimony(dataset, tree = start,
                      constraint = freeTaxaConstraint(),
                      maxReplicates = 4L, verbosity = 0L))
  expect_lt(attr(result, "score"), TreeLength(start, dataset))

  expect_true(all(vapply(result, SeparatesGroups, logical(1),
                         labels = labels, inGroup = c("a", "b"),
                         outGroup = c("c", "d"))))
})

test_that("the collapse pass keeps the enforced grouping visible", {
  # The one place the strict reading did return a wrong answer.  A constraint
  # is external evidence for a grouping, so ts_collapse_pool() protects the
  # branch that realises it from contraction — but it identified that branch by
  # matching a node's tip set to the `1` group EXACTLY.  With free taxa the
  # realising node is generally not that set, so nothing was protected and the
  # separating edge was contracted away: the returned tree broke the documented
  # constraint even though every tree the search visited satisfied it.
  labels <- letters[1:6]
  # Two characters support (a,e), two support (b,f); none supports (c,d), so
  # the branch that separates {a,b} from {c,d} is unsupported and collapses
  # unless it is protected.
  charDat <- StringToPhyDat(
    c("100010", "100010", "010001", "010001", "000000"), labels)
  at <- attributes(charDat)
  scoringConfig <- list(
    min_steps = integer(0), concavity = Inf, xpiwe = FALSE,
    xpiwe_r = 0.5, xpiwe_max_f = 5.0, obs_count = integer(0),
    infoAmounts = NULL
  )
  # {a,e,b,f} | {c,d} separates {a,b} from {c,d}; neither group is a clade.
  tree <- Preorder(RenumberTips(
    ape::read.tree(text = "(((a,e),(b,f)),(c,d));"), labels))
  cons <- phangorn::phyDat(
    matrix(c("1", "1", "0", "0", "?", "?"), nrow = 6,
           dimnames = list(labels, NULL)),
    type = "USER", levels = c("0", "1"))

  collapsed <- TreeSearch:::ts_collapse_pool(
    list(tree[["edge"]]), at$contrast,
    matrix(unlist(charDat, use.names = FALSE), nrow = 6, byrow = TRUE),
    at$weight, at$levels, scoringConfig, NULL, NULL,
    TreeSearch:::.PrepareConstraint(cons, charDat)[["consSplitMatrix"]])

  out <- structure(
    list(edge = collapsed$trees[[1]], tip.label = labels,
         Nnode = max(collapsed$trees[[1]]) - 6L),
    class = "phylo")
  expect_true(SeparatesGroups(Renumber(out), labels, c("a", "b"), c("c", "d")))
})

test_that(".PrepareConstraint codes free taxa as NA and drops vacuous rows", {
  dataset <- freeTaxaData()

  consArgs <- TreeSearch:::.PrepareConstraint(freeTaxaConstraint(), dataset)
  expect_equal(nrow(consArgs[["consSplitMatrix"]]), 1L)
  expect_equal(as.vector(consArgs[["consSplitMatrix"]]),
               c(1L, 1L, 0L, 0L, NA, NA, NA, NA))

  # A taxon on the tree but absent from the constraint is free too.
  partial <- phangorn::phyDat(
    matrix(c("1", "1", "0", "0"), nrow = 4,
           dimnames = list(letters[1:4], NULL)),
    type = "USER", levels = c("0", "1"))
  consArgs <- TreeSearch:::.PrepareConstraint(partial, dataset)
  expect_equal(as.vector(consArgs[["consSplitMatrix"]]),
               c(1L, 1L, 0L, 0L, NA, NA, NA, NA))

  # A group of fewer than two taxa is separated from the rest by every tree, so
  # such a character constrains nothing under the documented contract and must
  # not be enforced as a clade.  It is dropped, but not silently: coding only
  # `1` and `?` almost always means "group these taxa", which is not what it
  # says, and the alternative reading is the one that froze replicates.
  Inert <- function(...) {
    phangorn::phyDat(matrix(c(...), nrow = 8,
                            dimnames = list(letters[1:8], NULL)),
                     type = "USER", levels = c("0", "1"))
  }
  # No `0` group at all.
  expect_warning(
    dropped <- TreeSearch:::.PrepareConstraint(
      Inert("1", "1", "?", "?", "?", "?", "?", "?"), dataset),
    "trivial constraint")
  expect_equal(dropped, list())
  # A `0` group of one.  The two groups are interchangeable, so this must be
  # treated exactly like its mirror image below -- which the old
  # `1`-group-only test did not do.
  expect_warning(
    TreeSearch:::.PrepareConstraint(
      Inert("1", "1", "0", "?", "?", "?", "?", "?"), dataset),
    "trivial constraint")
  expect_warning(
    TreeSearch:::.PrepareConstraint(
      Inert("0", "0", "1", "?", "?", "?", "?", "?"), dataset),
    "trivial constraint")
  # Two and two: kept, and kept silently.
  expect_silent(TreeSearch:::.PrepareConstraint(
    Inert("1", "1", "0", "0", "?", "?", "?", "?"), dataset))

  # The loudest case of all, and the one the group-size test never sees: a
  # constraint with a single state, which is what MatrixToPhyDat() returns for
  # the `c(a = 1, b = 1, c = 1)` "make these a clade" idiom.  It reaches the
  # `nConsStates < 2` early return, so it must warn there.
  expect_warning(
    TreeSearch:::.PrepareConstraint(
      TreeTools::MatrixToPhyDat(c(a = "1", b = "1", c = "1")), dataset),
    "empty constraint")
})

test_that("the Wagner build places free taxa freely", {
  # wagner_tree_displays_constraint() and wagner_collect_active_splits() are a
  # second, independent implementation of the same reading.  AdditionTree() is
  # the path with no post-hoc retry to fall back on (has_posthoc is set only at
  # the search entry), so a Wagner build that read the constraint strictly
  # would warn here — and, before the fix, was forced to place every `?` taxon
  # outside the constrained group.
  dataset <- freeTaxaData()
  labels <- names(dataset)
  cons <- freeTaxaConstraint()

  for (seed in 1:8) {
    set.seed(seed)
    tree <- expect_silent(AdditionTree(dataset, constraint = cons))
    expect_true(SeparatesGroups(tree, labels, c("a", "b"), c("c", "d")),
                info = paste("seed", seed))
  }

  # A free taxon is genuinely free.  `wagner_collect_active_splits()` used to
  # read "outside the split" as ~split_tips, which put every `?` taxon in the
  # apart group and forced it out of the constrained clade: the tightest node
  # covering {a,b} and avoiding {c,d} was EXACTLY {a,b} in 25 of 25 seeds.  It
  # now holds at least one free taxon in all 25.  Asserting only that {a,b} and
  # {c,d} end up separated would not detect this -- an exact {a,b} clade
  # separates them too.
  tightest <- vapply(1:12, function(seed) {
    set.seed(seed)
    splits <- as.logical(as.Splits(AdditionTree(dataset, constraint = cons),
                                   tipLabels = labels))
    isOne <- labels %in% c("a", "b")
    isZero <- labels %in% c("c", "d")
    sizes <- c(
      rowSums(splits)[apply(splits, 1, function(r) {
        all(r[isOne]) && !any(r[isZero])
      })],
      (length(labels) - rowSums(splits))[apply(splits, 1, function(r) {
        !any(r[isOne]) && all(r[isZero])
      })]
    )
    if (length(sizes)) min(sizes) else NA_integer_
  }, numeric(1))
  # Compliant in every seed (no NA), and never pinned to the bare `1` group.
  expect_false(anyNA(tightest))
  expect_true(all(tightest > 2))
})
