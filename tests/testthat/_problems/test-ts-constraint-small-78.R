# Extracted from test-ts-constraint-small.R:78

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
