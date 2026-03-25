## Tests for constraint enforcement on small trees (T-039 regression)

library("TreeTools")

# Small 5-tip dataset for constraint tests
make_ds5 <- function() {
  phangorn::phyDat(
    matrix(c("0", "0", "0", "1", "1",
             "0", "1", "0", "1", "0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1")
  )
}

# Check that all constraint splits are displayed by the tree.
# Avoids %in%.Splits which has S3 dispatch issues in testthat's cloned
# namespace (test_check / R CMD check).
check_constraint <- function(tree, constraint) {
  tips <- sort(constraint$tip.label)
  tree_sp <- as.Splits(tree, tipLabels = tips)
  cons_sp <- as.Splits(constraint, tipLabels = tips)
  tm <- as.logical(tree_sp)
  cm <- as.logical(cons_sp)
  if (!is.matrix(tm)) tm <- matrix(tm, nrow = 1)
  if (!is.matrix(cm)) cm <- matrix(cm, nrow = 1)
  # Each constraint split must match some tree split (or its complement).
  # Use all()==/!= instead of identical() to avoid matrix vs vector mismatch.
  all(apply(cm, 1, function(c_row) {
    any(apply(tm, 1, function(t_row) {
      all(c_row == t_row) || all(c_row == !t_row)
    }))
  }))
}

test_that("T-039: fully resolving constraint on 5 tips doesn't crash", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(4217)
  expect_no_error(
    result <- MaximizeParsimony(ds5, constraint = cons,
                                maxReplicates = 1L, verbosity = 0L)
  )
  expect_s3_class(result, "multiPhylo")
  expect_equal(NTip(result[[1]]), 5L)
})

test_that("constraint satisfied on output trees (5 tips, 2 splits)", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  for (s in c(1, 7, 42, 99, 2718)) {
    set.seed(s)
    result <- MaximizeParsimony(ds5, constraint = cons,
                                maxReplicates = 2L, verbosity = 0L)
    for (i in seq_along(result)) {
      expect_true(
        check_constraint(result[[i]], cons),
        info = paste("seed", s, "tree", i)
      )
    }
  }
})

test_that("single constraint split on 5 tips works", {
  ds5 <- make_ds5()
  cons1 <- ape::read.tree(text = "((t1,t2),t3,t4,t5);")

  set.seed(3901)
  result <- MaximizeParsimony(ds5, constraint = cons1,
                              maxReplicates = 1L, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(check_constraint(result[[1]], cons1))
})

test_that("fully resolving constraint on 6 tips works", {
  ds6 <- phangorn::phyDat(
    matrix(c("0", "0", "0", "1", "1", "1",
             "0", "1", "0", "1", "0", "1"),
           nrow = 6, dimnames = list(paste0("t", 1:6), NULL)),
    type = "USER", levels = c("0", "1")
  )
  cons6 <- ape::read.tree(text = "((t1,t2),(t3,(t4,(t5,t6))));")

  set.seed(5537)
  result <- MaximizeParsimony(ds6, constraint = cons6,
                              maxReplicates = 1L, verbosity = 0L)
  expect_true(check_constraint(result[[1]], cons6))
})

test_that("constraint on 5 tips with IW scoring works", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(6614)
  result <- MaximizeParsimony(ds5, constraint = cons, concavity = 10,
                              maxReplicates = 1L, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(check_constraint(result[[1]], cons))
})

test_that("multiple replicates with constraint on small tree", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(8442)
  result <- MaximizeParsimony(ds5, constraint = cons,
                              maxReplicates = 5L, targetHits = 3L,
                              verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons),
                info = paste("tree", i))
  }
})

test_that("different fully-resolving constraints on 5 tips", {
  ds5 <- make_ds5()

  constraints <- list(
    ape::read.tree(text = "((t1,t3),(t2,(t4,t5)));"),
    ape::read.tree(text = "((t1,t4),(t2,(t3,t5)));"),
    ape::read.tree(text = "(t1,(t2,(t3,(t4,t5))));")
  )

  for (ci in seq_along(constraints)) {
    set.seed(1000 + ci)
    result <- MaximizeParsimony(ds5, constraint = constraints[[ci]],
                                maxReplicates = 1L, verbosity = 0L)
    expect_true(
      check_constraint(result[[1]], constraints[[ci]]),
      info = paste("constraint", ci)
    )
  }
})

test_that("T-208: adaptiveStart with constraints respects constraint", {
  ds5 <- make_ds5()
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(7293)
  result <- MaximizeParsimony(ds5, constraint = cons,
                              maxReplicates = 8L, targetHits = 4L,
                              verbosity = 0L,
                              control = SearchControl(adaptiveStart = TRUE))
  expect_s3_class(result, "multiPhylo")
  for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons),
                info = paste("tree", i, "violates constraint"))
  }
})
