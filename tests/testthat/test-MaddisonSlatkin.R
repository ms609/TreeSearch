library("TreeTools", quietly = TRUE)

test_that("MaddisonSlatkin() recursion bottoms", {
  expect_equal(MaddisonSlatkin(0, c(1, 1)), log(0))
  
  expect_equal(MaddisonSlatkin(1, c(1, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(2, 0)), log(1))
  expect_equal(MaddisonSlatkin(1, c(1, 0, 0, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(0, 0, 0, 2)), log(1))
})


test_that("MaddisonSlatkin() is numerically correct", {
  library("TreeTools", quietly = TRUE)
  
  expect_equal(MaddisonSlatkin(0, c(2, 1)), log(0))
  expect_equal(MaddisonSlatkin(1, c(2, 1)), log(1))
  expect_equal(MaddisonSlatkin(2, c(2, 1)), log(0))
  
  expect_slatkin <- function(tokens) {
    ch <- rep(seq_along(tokens), tokens)
    nTaxa <- length(ch)
    phyChar <- StringToPhyDat(paste0(ch, collapse = ""))
    trees <- as.phylo(1:NUnrooted(nTaxa) - 1, nTaxa)
    counts <- vapply(trees, TreeLength, double(1), phyChar) |>
      tabulate()
    out <- vapply(seq_along(counts), MaddisonSlatkin, double(1), tabulate(ch)) |>
      exp() * length(trees)
    expect_equal(out, counts)
  }
  expect_slatkin(c(2, 2))
  expect_slatkin(c(2, 3))
  expect_slatkin(c(2, 4))
  expect_slatkin(c(2, 3, 0, 2))
  
  exp(MaddisonSlatkin(2, c(2, 2)))
  exp(MaddisonSlatkin(1, c(2, 2)))
  LogCarter1(1,2,2)
  LogCarter1(2,2,2)
  
  
  # Maddison & Slatkin's tests
  expect_equal(MaddisonSlatkin(1, c(8, 24)) + LnUnrooted(32),
               LogCarter1(1, 8, 24))
  expect_equal(MaddisonSlatkin(2, c(8, 24)) + LnUnrooted(32),
               LogCarter1(2, 8, 24))
  
  # And a less even one
  expect_equal(MaddisonSlatkin(3, c(7, 18)) + LnUnrooted(25),
               LogCarter1(3, 7, 18))
})

test_that("MaddisonSlatkin() handles 5 states", {
  expect_equal(MaddisonSlatkin(4, c(2, 2, 0, 2, 0, 0, 0, 2, rep(0, 7), 2)),
               -6.851185) # by observation, not calculation
})
