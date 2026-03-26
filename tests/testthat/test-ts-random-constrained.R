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
