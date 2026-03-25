# Extracted from test-rearrange.cpp.R:5

# prequel ----------------------------------------------------------------------
library("TreeTools")

# test -------------------------------------------------------------------------
tr <- Preorder(RootTree(TreeTools::BalancedTree(7), "t1"))
expect_equal(0, length(expect_warning(all_tbr(tr$edge, -1))))
