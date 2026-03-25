## T-214: Multi-split constraints on 10+ tip trees
## Regression test for TBR rerooting destroying constraint splits
## that were classified as UNCONSTRAINED during clip phase.

# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

library("TreeTools")

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

test_that("T-214: two constraint splits on 10 tips", {
  set.seed(7142)
  m <- matrix(sample(c("0", "1"), 10 * 8, replace = TRUE),
              nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  ds10 <- phangorn::phyDat(m, type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2,t3),(t4,t5),t6,t7,t8,t9,t10);")

  for (s in c(137L, 274L, 411L, 548L, 685L)) {
    set.seed(s)
    result <- MaximizeParsimony(ds10, constraint = cons,
                                maxReplicates = 2L, verbosity = 0L)
    for (i in seq_along(result)) {
      expect_true(
        check_constraint(result[[i]], cons),
        info = paste("seed", s, "tree", i)
      )
    }
  }
})

test_that("T-214: three constraint splits on 15 tips", {
  set.seed(3241)
  m <- matrix(sample(c("0", "1"), 15 * 10, replace = TRUE),
              nrow = 15, dimnames = list(paste0("t", 1:15), NULL))
  ds15 <- phangorn::phyDat(m, type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(
    text = "((t1,t2,t3),(t4,t5,t6),(t7,t8),t9,t10,t11,t12,t13,t14,t15);"
  )

  for (s in c(43L, 86L, 129L, 172L, 215L)) {
    set.seed(s)
    result <- MaximizeParsimony(ds15, constraint = cons,
                                maxReplicates = 2L, verbosity = 0L)
    for (i in seq_along(result)) {
      expect_true(
        check_constraint(result[[i]], cons),
        info = paste("seed", s, "tree", i)
      )
    }
  }
})

test_that("T-214: nested constraint splits on 12 tips", {
  set.seed(5513)
  m <- matrix(sample(c("0", "1"), 12 * 6, replace = TRUE),
              nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  ds12 <- phangorn::phyDat(m, type = "USER", levels = c("0", "1"))
  # Nested: (t1,t2) inside (t1,t2,t3,t4)
  cons <- ape::read.tree(
    text = "((t1,t2),(t3,t4),t5,t6,t7,t8,t9,t10,t11,t12);"
  )

  for (s in c(311L, 622L, 933L)) {
    set.seed(s)
    result <- MaximizeParsimony(ds12, constraint = cons,
                                maxReplicates = 3L, verbosity = 0L)
    for (i in seq_along(result)) {
      expect_true(
        check_constraint(result[[i]], cons),
        info = paste("seed", s, "tree", i)
      )
    }
  }
})

test_that("T-214: multi-split constraint with IW scoring", {
  set.seed(7142)
  m <- matrix(sample(c("0", "1"), 10 * 8, replace = TRUE),
              nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  ds10 <- phangorn::phyDat(m, type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2,t3),(t4,t5),t6,t7,t8,t9,t10);")

  set.seed(2718)
  result <- MaximizeParsimony(ds10, constraint = cons, concavity = 10,
                              maxReplicates = 2L, verbosity = 0L)
  for (i in seq_along(result)) {
    expect_true(
      check_constraint(result[[i]], cons),
      info = paste("IW tree", i)
    )
  }
})
