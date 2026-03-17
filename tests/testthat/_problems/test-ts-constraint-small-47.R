# Extracted from test-ts-constraint-small.R:47

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TreeSearch", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library("TreeTools")
make_ds5 <- function() {
  phangorn::phyDat(
    matrix(c("0", "0", "0", "1", "1",
             "0", "1", "0", "1", "0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1")
  )
}
check_constraint <- function(tree, constraint) {
  tree_sp <- as.Splits(tree)
  cons_sp <- as.Splits(constraint)
  all(cons_sp %in% tree_sp)
}

# test -------------------------------------------------------------------------
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
