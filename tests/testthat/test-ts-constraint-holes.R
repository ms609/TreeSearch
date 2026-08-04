# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

## Three holes through which a `constraint` stopped binding the trees the user
## is handed (T-402, T-324, T-403).
##
## T-402: a start tree supplied via `tree =` was never checked against the
## constraint.  Constrained rearrangement cannot repair such a tree — an
## unmapped split makes regraft_violates_constraint() reject every move — so
## the replicate froze on it and reported its unconstrained score, which then
## evicted the compliant trees other replicates found.
##
## T-324: the per-replicate pool capture had no constraint gate at all,
## asymmetrically to the fuse capture beside it, so any violating tree that
## reached it was handed straight back.
##
## T-403: the collapse pass protected only a node whose tip set was the "1"
## group EXACTLY.  Tips ambiguous for a constraint character are free to sit on
## either side, so the split is often realised by a wider node — left
## unprotected, and contracted away under the default `collapse = TRUE`.

library("TreeTools", quietly = TRUE)

taxa <- letters[1:8]

# Does `tr` display a split with all of `one` on one side and all of `zero` on
# the other?  Spelled out from the edge matrix rather than via `Splits`: `%in%`
# on a Splits object dispatches differently under test_check() than under
# load_all(), and this has to answer the same way in both.
ConstraintShown <- function(tr, one, zero) {
  tr <- Postorder(tr)
  edge <- tr[["edge"]]
  label <- tr[["tip.label"]]
  nOne <- integer(max(edge))
  nZero <- integer(max(edge))
  nOne[match(one, label)] <- 1L
  nZero[match(zero, label)] <- 1L
  for (i in seq_len(nrow(edge))) {
    nOne[edge[i, 1]] <- nOne[edge[i, 1]] + nOne[edge[i, 2]]
    nZero[edge[i, 1]] <- nZero[edge[i, 1]] + nZero[edge[i, 2]]
  }
  # Postorder lists every node before its parent, so the first node holding a
  # whole group is that group's MRCA.
  nodes <- c(edge[, 2], edge[nrow(edge), 1])
  mrcaOne <- nodes[nOne[nodes] == length(one)][1]
  mrcaZero <- nodes[nZero[nodes] == length(zero)][1]
  (!is.na(mrcaOne) && nZero[mrcaOne] == 0L) ||
    (!is.na(mrcaZero) && nOne[mrcaZero] == 0L)
}

AllShown <- function(trees, one, zero) {
  sum(vapply(trees, ConstraintShown, logical(1), one, zero))
}

# Two characters supporting (a, e) and two supporting (b, f): the unconstrained
# optimum groups a with e and b with f, which no tree holding {a, b} together
# can do.  Optimum 6 unconstrained, 10 under an {a, b} constraint.
abDataset <- local({
  m <- rbind(
    c(1, 0, 0, 0, 1, 0, 0, 0),
    c(1, 0, 0, 0, 1, 0, 0, 0),
    c(0, 1, 0, 0, 0, 1, 0, 0),
    c(0, 1, 0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 0, 0, 1, 1),
    c(0, 0, 0, 0, 0, 0, 1, 1)
  )
  colnames(m) <- taxa
  MatrixToPhyDat(t(m))
})

# {a, b} against every other taxon: no ambiguous tip, so the constraint the
# user states and the stricter one the search enforces internally coincide.
abConstraint <- MatrixToPhyDat(matrix(
  c(1, 1, 0, 0, 0, 0, 0, 0), ncol = 1, dimnames = list(taxa, NULL)
))

# Scores 6 -- better than any {a, b}-compliant tree -- and violates {a, b}.
abViolatingStart <- ape::read.tree(text = "(((a,e),(b,f)),((c,d),(g,h)));")


test_that("a violating `tree` cannot beat the constrained optimum (T-402)", {
  # Control: the same constraint from a cold start reaches 10 and complies.
  set.seed(1)
  cold <- MaximizeParsimony(abDataset, constraint = abConstraint,
                            maxReplicates = 4L, verbosity = 0L)
  expect_equal(as.numeric(attr(cold, "score")), 10)
  expect_equal(AllShown(cold, c("a", "b"), setdiff(taxa, c("a", "b"))),
               length(cold))

  # maxReplicates = 1: gating the pool capture alone would leave the pool
  # empty here, and MaximizeParsimony() would fall back to returning the
  # supplied start.  The start must be dealt with at the boundary.
  set.seed(1)
  expect_warning(
    one <- MaximizeParsimony(abDataset, tree = abViolatingStart,
                             constraint = abConstraint, maxReplicates = 1L,
                             verbosity = 0L),
    "do not satisfy `constraint`"
  )
  expect_equal(as.numeric(attr(one, "score")), 10)
  expect_equal(AllShown(one, c("a", "b"), setdiff(taxa, c("a", "b"))),
               length(one))

  # Several replicates: the violating tree's illegal score used to evict every
  # compliant tree the other replicates found.
  set.seed(1)
  expect_warning(
    many <- MaximizeParsimony(abDataset, tree = abViolatingStart,
                              constraint = abConstraint, maxReplicates = 8L,
                              verbosity = 0L),
    "do not satisfy `constraint`"
  )
  expect_equal(as.numeric(attr(many, "score")), 10)
  expect_equal(AllShown(many, c("a", "b"), setdiff(taxa, c("a", "b"))),
               length(many))
})


test_that("a violating tree never enters the pool (T-324)", {
  # The Wagner retry-exhaustion route that motivated T-324 is not constructible
  # on demand, so the shared downstream half -- the ungated pool capture -- is
  # driven through T-402's start instead: without the gate the replicate's
  # frozen, violating tree is captured verbatim.  `poolSuboptimal` keeps
  # non-best trees too, so a violating tree would be visible even if a better
  # compliant one existed.
  set.seed(2)
  expect_warning(
    result <- MaximizeParsimony(abDataset, tree = abViolatingStart,
                                constraint = abConstraint, maxReplicates = 3L,
                                verbosity = 0L, collapse = FALSE,
                                poolSuboptimal = 4),
    "do not satisfy `constraint`"
  )
  expect_equal(AllShown(result, c("a", "b"), setdiff(taxa, c("a", "b"))),
               length(result))
  expect_gte(as.numeric(attr(result, "score")), 10)
})


test_that("collapse keeps the constraint visible (T-403)", {
  # Only (a, e) and (b, f) are supported, so the branch that separates
  # {a, b} from {c, d} is unsupported and collapses -- taking the constraint
  # with it.  The node realising the split is {a, e, b, f}, not the "1" group
  # {a, b}, which is why an exact-match protection missed it.
  m <- rbind(
    c(1, 0, 0, 0, 1, 0, 0, 0),
    c(1, 0, 0, 0, 1, 0, 0, 0),
    c(0, 1, 0, 0, 0, 1, 0, 0),
    c(0, 1, 0, 0, 0, 1, 0, 0)
  )
  colnames(m) <- taxa
  dataset <- MatrixToPhyDat(t(m))

  # e--h ambiguous: the constraint asks only that {a, b} be separated from
  # {c, d}, which the start below already does.
  constraint <- MatrixToPhyDat(matrix(
    c("1", "1", "0", "0", "?", "?", "?", "?"),
    ncol = 1, dimnames = list(taxa, NULL)
  ))
  start <- ape::read.tree(text = "(((a,e),(b,f)),(c,(d,(g,h))));")

  set.seed(1)
  collapsed <- MaximizeParsimony(dataset, tree = start,
                                 constraint = constraint, maxReplicates = 2L,
                                 verbosity = 0L)
  expect_equal(AllShown(collapsed, c("a", "b"), c("c", "d")),
               length(collapsed))

  # Built-in control: without collapsing, the split was never at risk.
  set.seed(1)
  resolved <- MaximizeParsimony(dataset, tree = start,
                                constraint = constraint, maxReplicates = 2L,
                                verbosity = 0L, collapse = FALSE)
  expect_equal(AllShown(resolved, c("a", "b"), c("c", "d")),
               length(resolved))
  expect_equal(as.numeric(attr(collapsed, "score")),
               as.numeric(attr(resolved, "score")))
})
