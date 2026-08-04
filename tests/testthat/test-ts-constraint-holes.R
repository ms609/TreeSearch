# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

## A `constraint` must bind the trees the caller is handed, at each of the three
## boundaries where it can be lost: the starting tree, the pool capture, and the
## final collapse (T-402, T-324, T-403).
##
## Each test asserts COMPLIANCE of the returned trees, not the score alone.  A
## constraint-violating tree is drawn from a wider set of topologies than a legal
## one, so it scores better; a score assertion alone would pass on exactly the
## tree that breaks the contract.

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

  # Several replicates: an illegal score is better than any legal one, so it
  # evicts every compliant tree the other replicates find.  One bad start must
  # not cost the whole search.
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


test_that("no tree in the returned pool breaks the constraint (T-324)", {
  # What this asserts is the outcome -- every tree handed back complies -- over
  # the whole pool, not just the best-score trees: `poolSuboptimal` retains the
  # near-misses, which is where an ungated capture shows up.  It does NOT prove
  # the capture gate itself fires; the route that motivated T-324 is Wagner
  # retry-exhaustion, whose reachability is unconfirmed and which cannot be
  # forced from R.  Treat this as a contract test, not a gate test.
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

  # The parallel driver has its own copy of the capture, on a per-thread
  # constraint and pool; two threads is the project's per-agent core limit.
  set.seed(2)
  expect_warning(
    parallel <- MaximizeParsimony(abDataset, tree = abViolatingStart,
                                  constraint = abConstraint, maxReplicates = 4L,
                                  nThreads = 2L, verbosity = 0L),
    "do not satisfy `constraint`"
  )
  expect_equal(AllShown(parallel, c("a", "b"), setdiff(taxa, c("a", "b"))),
               length(parallel))
  expect_equal(as.numeric(attr(parallel, "score")), 10)
})


test_that("every flat kernel takes .PrepareConstraint()'s output", {
  # The flat `ts_*` kernels declare their constraint arguments as formals, so a
  # field .PrepareConstraint() adds for the list-config entry points is an
  # unused-argument error at any site that splats the whole list into one.
  # Assert the filter covers every formal each kernel actually declares, and
  # exercise the entry points that splat -- `Resample(nReplicates > 1)` had no
  # constrained coverage at all, so an unfiltered splat there stayed green.
  kernels <- list(TreeSearch:::ts_wagner_tree,
                  TreeSearch:::ts_random_wagner_tree,
                  TreeSearch:::ts_resample_search,
                  TreeSearch:::ts_parallel_resample,
                  TreeSearch:::ts_successive_approx)
  filtered <- names(TreeSearch:::.KernelConstraintArgs(
    TreeSearch:::.PrepareConstraint(abConstraint, abDataset)
  ))
  for (k in kernels) {
    expect_true(all(filtered %in% names(formals(k))))
  }

  set.seed(4)
  expect_s3_class(
    Resample(abDataset, constraint = abConstraint, nReplicates = 2L,
             maxReplicates = 2L),
    "multiPhylo"
  )
  set.seed(4)
  expect_s3_class(AdditionTree(abDataset, constraint = abConstraint), "phylo")
})


test_that("collapse keeps the constraint visible (T-403)", {
  # Only (a, e) and (b, f) are supported, so the branch that separates
  # {a, b} from {c, d} is unsupported and collapses -- taking the constraint
  # with it.  The node realising the split is {a, e, b, f}, not the "1" group
  # {a, b}, so protection keyed on an exact match with the "1" group does not
  # reach it.
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
