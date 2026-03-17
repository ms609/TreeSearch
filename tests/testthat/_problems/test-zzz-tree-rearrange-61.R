# Extracted from test-zzz-tree-rearrange.R:61

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
testTree <- Preorder(root(BalancedTree(7), 1, resolve.root = TRUE))
edge <- testTree[["edge"]]
expect_equal(spr(edge, 66), cSPR(testTree, 66)$edge)
