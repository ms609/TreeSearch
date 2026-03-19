# Extracted from test-Morphy.R:7

# prequel ----------------------------------------------------------------------
library("TreeTools", quietly = TRUE)

# test -------------------------------------------------------------------------
skip_if(interactive())
dataset <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 3, f = 3))
expect_warning(PrepareDataProfile(dataset),
                 "Can handle max. 2 informative tokens")
