# Extracted from test-zzz-tree-rearrange.R:38

# prequel ----------------------------------------------------------------------
library("TreeTools")
context("Tree rearrangements")
tree5a <- read.tree(text = '(a, (b, (c, (d, e))));')
tree5b <- read.tree(text = '((a, b), (c, (d, e)));')
tree6  <- Preorder(read.tree(text = "((a, (b, (c, d))), (e, f));"))
tree6b <- Preorder(read.tree(text = "((a, (b, c)), (d, (e, f)));"))
tree8  <- read.tree(text = "(((a, (b, (c, d))), (e, f)), (g, h));")
tree11 <- read.tree(text = "((((a, b), (c, d)), e), ((f, (g, (h, i))), (j, k)));")
attr(tree5a, 'order') <- attr(tree5b, 'order') <- attr(tree8, 'order') <- attr(tree11, 'order') <- 'preorder'

# test -------------------------------------------------------------------------
trComb <- read.tree(text = "(((((1,2),3),4),5),6);")
edge <- trComb$edge
Test <- function (e, r, e1, e2) {
    edge1 <- edge
    edge1[c(e1, e2), 2] <- edge1[c(e2, e1), 2]
    edge1 <- do.call(cbind, RenumberEdges(edge1[, 1], edge1[, 2]))
    expect_equal(edge1, nni(trComb$edge, e, r))
  }
