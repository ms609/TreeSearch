# Extracted from test-rearrange.cpp.R:69

# prequel ----------------------------------------------------------------------
library("TreeTools")

# test -------------------------------------------------------------------------
t2 <- as.phylo(518, 7)
expect_equal(8, length(all_spr(t2$edge, 2)))
