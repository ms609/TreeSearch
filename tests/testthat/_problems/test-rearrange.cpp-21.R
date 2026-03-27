# Extracted from test-rearrange.cpp.R:21

# prequel ----------------------------------------------------------------------
library("TreeTools")

# test -------------------------------------------------------------------------
tr <- Preorder(RootTree(TreeTools::BalancedTree(7), "t1"))
expect_equal(8, length(x <- all_tbr(tr$edge, 12)))
