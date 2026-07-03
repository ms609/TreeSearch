## T-212: Test RANDOM_TREE strategy with constraints.
## Exercises random_constrained_tree() via ts_driven_search: builds the
## constraint backbone then randomly resolves polytomies. Tests verify
## constraint satisfaction across serial, parallel, adaptive-start, IW,
## and single-split scenarios.

skip_on_cran()
library("TreeTools")

make_ds5 <- function() {
  phangorn::phyDat(
    matrix(c("0", "0", "0", "1", "1",
             "0", "1", "0", "1", "0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1")
  )
}

# Constraint checker (same as test-ts-constraint-small.R)
check_constraint <- function(tree, constraint) {
  tips <- sort(constraint$tip.label)
  tree_sp <- as.Splits(tree, tipLabels = tips)
  cons_sp <- as.Splits(constraint, tipLabels = tips)
  tm <- as.logical(tree_sp)
  cm <- as.logical(cons_sp)
  if (!is.matrix(tm)) tm <- matrix(tm, nrow = 1)
  if (!is.matrix(cm)) cm <- matrix(cm, nrow = 1)
  all(apply(cm, 1, function(c_row) {
    any(apply(tm, 1, function(t_row) {
      all(c_row == t_row) || all(c_row == !t_row)
    }))
  }))
}

test_that("RANDOM_TREE strategy with constraint (serial, 5 tips)", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  for (s in c(3142L, 5263L, 7384L)) {
    set.seed(s)
    result <- MaximizeParsimony(
      ds5, constraint = cons,
      maxReplicates = 4L, verbosity = 0L,
      control = SearchControl(wagnerBias = 3L)
    )
    expect_s3_class(result, "multiPhylo")
    for (i in seq_along(result)) {
      expect_true(
        check_constraint(result[[i]], cons),
        info = paste("seed", s, "tree", i)
      )
    }
  }
})

test_that("RANDOM_TREE score is valid", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(8461)
  result <- MaximizeParsimony(
    ds5, constraint = cons,
    maxReplicates = 2L, verbosity = 0L,
    control = SearchControl(wagnerBias = 3L)
  )
  for (i in seq_along(result)) {
    score <- TreeLength(result[[i]], ds5)
    expect_true(score > 0, info = paste("tree", i, "score > 0"))
    expect_true(is.finite(score), info = paste("tree", i, "score finite"))
  }
})

test_that("adaptiveStart round-robin with constraints (serial)", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  # 8 replicates → bandit samples RANDOM_TREE for some
  set.seed(9527)
  result <- MaximizeParsimony(
    ds5, constraint = cons,
    maxReplicates = 8L, targetHits = 4L,
    verbosity = 0L,
    control = SearchControl(adaptiveStart = TRUE)
  )
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(
      check_constraint(result[[i]], cons),
      info = paste("adaptive serial, tree", i)
    )
  }
})

test_that("adaptiveStart round-robin with constraints (parallel, nThreads=2)", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  # 8 replicates with 2 threads → parallel round-robin.
  # RANDOM_TREE assigned at rep 3, 7 (r % 4 == 3).
  set.seed(1653)
  result <- MaximizeParsimony(
    ds5, constraint = cons,
    maxReplicates = 8L, targetHits = 4L,
    verbosity = 0L, nThreads = 2L,
    control = SearchControl(adaptiveStart = TRUE)
  )
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(
      check_constraint(result[[i]], cons),
      info = paste("adaptive parallel, tree", i)
    )
  }
})

test_that("RANDOM_TREE with single constraint split", {
  ds5 <- make_ds5()
  cons1 <- ape::read.tree(text = "((t1,t2),t3,t4,t5);")

  set.seed(6274)
  result <- MaximizeParsimony(
    ds5, constraint = cons1,
    maxReplicates = 4L, verbosity = 0L,
    control = SearchControl(wagnerBias = 3L)
  )
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(
      check_constraint(result[[i]], cons1),
      info = paste("single split, tree", i)
    )
  }
})

test_that("RANDOM_TREE with IW scoring + constraints", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(5018)
  result <- MaximizeParsimony(
    ds5, constraint = cons, concavity = 10,
    maxReplicates = 4L, verbosity = 0L,
    control = SearchControl(wagnerBias = 3L)
  )
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(
      check_constraint(result[[i]], cons),
      info = paste("IW + constraint, tree", i)
    )
  }
})

## T-329: an impossible (non-laminar) constraint must error cleanly rather
## than reach random_constrained_tree(), which cannot represent splits that
## no tree can display simultaneously.
test_that("impossible (non-laminar) constraint errors cleanly", {
  ds5 <- make_ds5()

  # Split A = {t1,t2,t3}, split B = {t3,t4,t5}: overlap {t3} is neither
  # empty nor equal to either split, so A and B cannot both appear on a
  # tree.
  consMat <- matrix(c("1", "1", "1", "0", "0",
                       "0", "0", "1", "1", "1"),
                     nrow = 5, dimnames = list(paste0("t", 1:5), NULL))
  cons <- phangorn::phyDat(consMat, type = "USER", levels = c("0", "1"))

  expect_error(
    TreeSearch:::.PrepareConstraint(cons, ds5),
    "impossible"
  )
  expect_error(
    MaximizeParsimony(ds5, constraint = cons,
                       maxReplicates = 1L, verbosity = 0L),
    "impossible"
  )
})

## T-329 (defensive fix): random_constrained_tree() must never inject a
## collapsed split's phantom node (-1) into its parent's item list.
##
## A split "collapses" (split_root == -1) when every tip it owns is stolen
## by tighter, non-laminar splits that are not its formal children.  If that
## collapsed split still has a strict-superset parent, the parent's child
## collection loop (src/ts_wagner.cpp, "Child splits whose parent is this
## split") must skip it rather than push -1.
##
## This constructs that scenario directly via the internal flat
## ts_driven_search() interface, bypassing .PrepareConstraint()'s R-level
## laminarity check (T-329's primary fix) so the C++ guard itself is
## exercised:
##   D1 = {t1,t6}, D2 = {t2,t7}, D3 = {t3,t8}   (tight, non-laminar splits)
##   C  = {t1,t2,t3}                            (loses every tip to D1-D3)
##   S  = {t1,t2,t3,t6,t7,t8}                   (strict superset of C, D1-D3)
## C collapses (split_root == -1) but is still S's formal child, so S's
## build loop must skip it.  wagnerBias = 3L forces the RANDOM_TREE
## strategy so random_constrained_tree() is actually invoked.
test_that("collapsed constraint split with superset parent does not corrupt tree", {
  tipNames <- c("t9", "t1", "t2", "t3", "t6", "t7", "t8")
  ds <- phangorn::phyDat(
    matrix(c("0", "1", "0", "1", "0", "1", "0",
             "1", "0", "1", "0", "1", "0", "1"),
           nrow = 7, dimnames = list(tipNames, NULL)),
    type = "USER", levels = c("0", "1")
  )
  at <- attributes(ds)

  consSplitMatrix <- matrix(c(
    0, 1, 0, 0, 1, 0, 0,   # D1 = {t1, t6}
    0, 0, 1, 0, 0, 1, 0,   # D2 = {t2, t7}
    0, 0, 0, 1, 0, 0, 1,   # D3 = {t3, t8}
    0, 1, 1, 1, 0, 0, 0,   # C  = {t1, t2, t3}  (collapses)
    0, 1, 1, 1, 1, 1, 1    # S  = {t1, t2, t3, t6, t7, t8}
  ), nrow = 5, byrow = TRUE)

  for (s in c(1L, 2L, 3L)) {
    set.seed(s)
    res <- TreeSearch:::ts_driven_search(
      contrast = at$contrast,
      tip_data = matrix(unlist(ds, use.names = FALSE), nrow = 7,
                         byrow = TRUE),
      weight = at$weight,
      levels = at$levels,
      maxReplicates = 1L, targetHits = 1L, verbosity = 0L,
      wagnerBias = 3L,
      consSplitMatrix = consSplitMatrix
    )
    edge <- res$trees[[1]]
    n_tip <- 7L
    n_node <- 2L * n_tip - 1L

    expect_equal(nrow(edge), 2L * n_tip - 2L, info = paste("seed", s))
    expect_true(all(edge >= 1L & edge <= n_node), info = paste("seed", s))
    # Every node 1..n_node bar the root appears as a child exactly once.
    expect_equal(sort(unique(edge[, 2])), setdiff(seq_len(n_node),
                                                    n_tip + 1L),
                 info = paste("seed", s))
  }
})

test_that("parallel RANDOM_TREE with multiple seeds (nThreads=2)", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  for (s in c(2241L, 4362L)) {
    set.seed(s)
    result <- MaximizeParsimony(
      ds5, constraint = cons,
      maxReplicates = 8L, targetHits = 4L,
      verbosity = 0L, nThreads = 2L,
      control = SearchControl(adaptiveStart = TRUE)
    )
    expect_s3_class(result, "multiPhylo")
    for (i in seq_along(result)) {
      expect_true(
        check_constraint(result[[i]], cons),
        info = paste("parallel seed", s, "tree", i)
      )
    }
  }
})
