## agent-issues/TreeSearch#54 (follow-up): random_constrained_tree() must
## sample the trees the constraint contract ALLOWS, not the corner of them in
## which every free taxon sits outside every constrained group.
##
## The contract (`?MaximizeParsimony`, `@param constraint`) is that a tree
## complies when some edge separates the taxa coded 1 from those coded 0,
## `?`-coded taxa falling on either side.  The generator built its backbone from
## the "together" group alone and made every free tip a root-level item, so the
## group always came out as an EXACT clade and no free tip ever started inside
## it.  Compliant, but only a fraction of the compliant trees were reachable.

skip_on_cran()
library("TreeTools")

## Edge matrix (as returned by the C++ generator) -> phylo.
rctPhylo <- function(edge, tips) {
  structure(list(edge = edge, tip.label = tips,
                 Nnode = length(tips) - 1L),
            class = "phylo")
}

## Split membership matrix, columns in `tips` order.
rctSplits <- function(tree, tips) {
  sp <- as.Splits(tree, tipLabels = tips)
  m <- as.logical(sp)
  if (!is.matrix(m)) m <- matrix(m, nrow = 1)
  colnames(m) <- attr(sp, "tip.label")
  m[, tips, drop = FALSE]
}

## Does some edge put all of `together` on one side and all of `apart` on the
## other?  This is the documented contract, stated without reference to which
## group the machinery happens to canonicalise as "inside".
rctSeparates <- function(tree, tips, together, apart) {
  m <- rctSplits(tree, tips)
  any(apply(m, 1, function(r) {
    all(r[together] == r[together][[1]]) &&
      all(r[apart] == r[apart][[1]]) &&
      r[together][[1]] != r[apart][[1]]
  }))
}

## Size of the smallest clade that holds every tip of `together` and none of
## `apart`; NA if no edge separates them.  Equals length(together) exactly when
## the group is an exact clade, i.e. when no free tip sits inside it.
rctTightest <- function(tree, tips, together, apart) {
  m <- rctSplits(tree, tips)
  sizes <- apply(m, 1, function(r) {
    for (side in list(r, !r)) {
      if (all(side[together]) && !any(side[apart])) return(sum(side))
    }
    NA_integer_
  })
  sizes <- sizes[!is.na(sizes)]
  if (length(sizes)) min(sizes) else NA_integer_
}

rctDataset <- function(tips) {
  n <- length(tips)
  phangorn::phyDat(
    matrix(c(rep_len(c("0", "1"), n), rep_len(c("1", "0"), n)),
           nrow = n, dimnames = list(tips, NULL)),
    type = "USER", levels = c("0", "1")
  )
}

rctDraw <- function(tsd, splitMatrix, tips) {
  rctPhylo(
    TreeSearch:::ts_random_constrained_tree(
      tsd$contrast, tsd$tip_data, tsd$weight, tsd$levels,
      consSplitMatrix = splitMatrix),
    tips)
}

## Exhaustive: 6 taxa is small enough to enumerate every unrooted binary tree
## and say exactly which ones the constraint allows, so "does it sample at
## random from the legal set" has a yes/no answer rather than a distributional
## one.  35 of the 105 trees separate {c,d} from {a,b}.  Before the fix the
## generator could return only 15 of them; the other 20 were unreachable at
## every seed.
test_that("random_constrained_tree samples every legal topology", {
  tips <- letters[1:6]
  ds <- rctDataset(tips)
  tsd <- make_ts_data(ds)
  # {c,d} together, {a,b} apart, {e,f} free.  Coded 1/0/NA exactly as
  # .PrepareConstraint() writes it.
  splitMatrix <- matrix(c(0L, 0L, 1L, 1L, NA_integer_, NA_integer_), nrow = 1)

  treeNo <- function(tree) as.character(as.numeric(as.TreeNumber(tree)))
  everyTree <- lapply(seq_len(105) - 1,
                      function(i) as.phylo(i, nTip = 6, tipLabels = tips))
  legal <- vapply(everyTree, rctSeparates, logical(1),
                  tips = tips, together = c("c", "d"), apart = c("a", "b"))
  expect_equal(sum(legal), 35L)
  legalNos <- vapply(everyTree[legal], treeNo, character(1))

  seen <- character(2000)
  for (s in seq_along(seen)) {
    set.seed(s)
    seen[[s]] <- treeNo(rctDraw(tsd, splitMatrix, tips))
  }

  # Sound: never a tree the constraint forbids.
  expect_equal(setdiff(seen, legalNos), character(0))
  # Complete: every tree the constraint permits is reachable.  This is what
  # fails pre-fix -- 20 of the 35 never appear.
  expect_equal(sort(setdiff(legalNos, seen)), character(0))
  # ...and reachable at a comparable rate, not merely grazed.  Expectation is
  # 2000 / 35 = 57 draws each.
  hits <- as.vector(table(factor(seen, levels = legalNos)))
  expect_gt(min(hits), 20)
  expect_lt(max(hits), 120)
})

## The user-facing route in: `?`-coded taxa in a constraint phyDat, through
## .PrepareConstraint(), with the group the machinery canonicalises as "inside"
## chosen by tip 0's coding rather than by the test.
test_that("`?` taxa start inside a constrained group as well as outside", {
  tips <- letters[1:8]
  ds <- rctDataset(tips)
  cons <- phangorn::phyDat(
    matrix(c("1", "1", "0", "0", "?", "?", "?", "?"),
           ncol = 1, dimnames = list(tips, NULL)),
    type = "USER", levels = c("0", "1")
  )
  splitMatrix <- TreeSearch:::.PrepareConstraint(cons, ds)$consSplitMatrix
  tsd <- make_ts_data(ds)

  # build_constraint() swaps the groups so that tip 0 (a) is never in the
  # "inside" mask, so {c,d} is the group this builds as a clade.  That is the
  # one a free tip could never join.
  tightest <- integer(50)
  for (s in seq_along(tightest)) {
    set.seed(s)
    tree <- rctDraw(tsd, splitMatrix, tips)
    expect_true(rctSeparates(tree, tips, c("a", "b"), c("c", "d")),
                info = paste("seed", s))
    tightest[[s]] <- rctTightest(tree, tips, c("c", "d"), c("a", "b"))
  }
  # Pre-fix this is 2 at every seed: {c,d} is always exactly a clade.
  expect_gt(max(tightest), 2L)
  expect_equal(min(tightest), 2L)  # and still sometimes exactly a clade
})

## Guard against over-loosening: a constraint that names every taxon has no
## free tips, so the generator must behave exactly as it always did.
test_that("a constraint with no free taxa still builds exact clades", {
  tips <- letters[1:6]
  ds <- rctDataset(tips)
  tsd <- make_ts_data(ds)
  splitMatrix <- matrix(c(0L, 0L, 1L, 1L, 0L, 0L), nrow = 1)

  for (s in 1:25) {
    set.seed(s)
    tree <- rctDraw(tsd, splitMatrix, tips)
    expect_equal(
      rctTightest(tree, tips, c("c", "d"), c("a", "b", "e", "f")), 2L,
      info = paste("seed", s)
    )
  }
})

## Structural validity with free tips present.  Scattering them consumes
## internal node indices that the backbone did not, so the node budget
## (n_tip - 1 internal nodes, none allocated twice, none left dangling) is worth
## pinning: an over-run would corrupt the tree rather than fail loudly.
## Two splits, and a free tip 0 -- the tip whose position build_constraint()
## canonicalises the split masks around.
test_that("scattered free tips leave a well-formed tree", {
  tips <- paste0("t", 1:9)
  ds <- rctDataset(tips)
  tsd <- make_ts_data(ds)
  nTip <- 9L
  nNode <- 2L * nTip - 1L
  # t1 (tip 0) free in both splits; {t2,t3} vs {t4,t5}; {t6,t7} vs {t8,t9}.
  na <- NA_integer_
  splitMatrix <- matrix(c(
    na, 1L, 1L, 0L, 0L, na, na, na, na,
    na, na, na, na, na, 1L, 1L, 0L, 0L
  ), nrow = 2, byrow = TRUE)

  for (s in 1:25) {
    set.seed(s)
    edge <- TreeSearch:::ts_random_constrained_tree(
      tsd$contrast, tsd$tip_data, tsd$weight, tsd$levels,
      consSplitMatrix = splitMatrix)
    expect_equal(nrow(edge), 2L * nTip - 2L, info = paste("seed", s))
    expect_true(all(edge >= 1L & edge <= nNode), info = paste("seed", s))
    # Every node but the root is somebody's child, exactly once.
    expect_equal(sort(edge[, 2]), setdiff(seq_len(nNode), nTip + 1L),
                 info = paste("seed", s))
    tree <- rctPhylo(edge, tips)
    expect_true(rctSeparates(tree, tips, c("t2", "t3"), c("t4", "t5")),
                info = paste("seed", s))
    expect_true(rctSeparates(tree, tips, c("t6", "t7"), c("t8", "t9")),
                info = paste("seed", s))
  }
})
