# Extracted from test-rearrange.cpp.R:21

# prequel ----------------------------------------------------------------------
library("TreeTools")

# test -------------------------------------------------------------------------
tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
expect_equal(8, length(x <- all_tbr(tr$edge, 12)))
