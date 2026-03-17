# Extracted from test-Morphy.R:193

# prequel ----------------------------------------------------------------------
library("TreeTools", quietly = TRUE)

# test -------------------------------------------------------------------------
x <- structure(
    array(c(
      rep(1L, 8),
      rep(2L, 8),
      rep(3L, 8),
      rep(2L, 8),
      rep(1L, 8)
      ),
      dim = c(4, 2, 5)),
    firstHit = c(start = 5, test = 0, end = 0)
  )
y <- array(c(rep(1L, 8),
               rep(4L, 8),
               rep(1L, 8),
               rep(4L, 8),
               rep(1L, 8)),
          dim = c(4, 2, 5)
          )
expect_warning(.CombineResults(x, y, stage = "test"))
