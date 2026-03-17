# Extracted from test-Concordance.R:123

# prequel ----------------------------------------------------------------------
library("TreeTools", quietly = TRUE)

# test -------------------------------------------------------------------------
expect_equal(.Rezero(seq(0, 1, by = 0.1), 0.1), -1:9 / 9)
