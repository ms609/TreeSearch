# Extracted from test-ts-constraint-small.R:100

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
set.seed(8442)
result <- MaximizeParsimony(ds5, constraint = cons,
                              maxReplicates = 5L, targetHits = 3L,
                              verbosity = 0L)
expect_s3_class(result, "multiPhylo")
for (i in seq_along(result)) {
    expect_true(check_constraint(result[[i]], cons),
                info = paste("tree", i))
  }
