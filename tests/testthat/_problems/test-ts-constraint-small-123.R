# Extracted from test-ts-constraint-small.R:123

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
