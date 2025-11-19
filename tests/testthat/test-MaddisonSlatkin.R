test_that(".LogR() matches empirical observation", {
  # R is the probability that in a randomly selected tree of _n_ taxa, the
  # smaller of the two basal subclades will have _m_ taxa.
  library("TreeTools", quietly = TRUE)
  
  n <- 6
  # All trees are rooted on t1, which we'll therefore use as the root
  trees <- as.phylo(1:NUnrooted(n + 1) - 1, n + 1)
  counts <- vapply(trees, function(tree) {
    ed <- tree$edge
    min(CladeSizes(tree, nodes = ed[ed[, 1] == n + 3, 2]))
    }, double(1)) |>
    table() |>
    as.numeric()
  expect_equal(exp(c(.LogR(1, 6), .LogR(2, 6), .LogR(3, 6))) * NRooted(n),
               counts)
})

test_that(".LogD() succeeds", {
  # D is the probability that, in a randomly selected tree on `leaves`, the
  # smaller subclade of taxa will receive taxa with labels `drawn`
  
})



test_that("MaddisonSlatkin() recursion bottoms", {
  expect_equal(MaddisonSlatkin(0, c(1, 1)), log(0))
  
  expect_equal(MaddisonSlatkin(1, c(1, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(2, 0)), log(1))
  expect_equal(MaddisonSlatkin(1, c(1, 0, 0, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(0, 0, 0, 2)), log(1))
  
  
test_that("MaddisonSlatkin is numerically correct", {
  expect_equal(MaddisonSlatkin(0, c(2, 1)), log(0))
  expect_equal(MaddisonSlatkin(1, c(2, 1)), TreeTools::LnRooted(3))
  expect_equal(MaddisonSlatkin(2, c(2, 1)), log(0))
})
