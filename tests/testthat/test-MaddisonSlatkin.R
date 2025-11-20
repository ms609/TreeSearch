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
  expect_equal(exp(c(.LogR(1, n), .LogR(2, n), .LogR(3, n))) * NRooted(n),
               counts)
  
  n <- 7
  trees <- as.phylo(1:NUnrooted(n + 1) - 1, n + 1)
  counts <- vapply(trees, function(tree) {
    ed <- tree$edge
    min(CladeSizes(tree, nodes = ed[ed[, 1] == n + 3, 2]))
    }, double(1)) |>
    table() |>
    as.numeric()
  expect_equal(exp(c(.LogR(1, n), .LogR(2, n), .LogR(3, n))) * NRooted(n),
               counts)
  
  # Edge cases
  n <- 2
  trees <- c(as.phylo(1:NUnrooted(n + 1) - 1, n + 1))
  counts <- vapply(trees, function(tree) {
    ed <- tree$edge
    min(CladeSizes(tree, nodes = ed[ed[, 1] == n + 3, 2]))
  }, double(1)) |>
    table() |>
    as.numeric()
  expect_equal(exp(c(.LogR(1, n))) * NRooted(n), counts)
  
  n <- 3
  trees <- as.phylo(1:NUnrooted(n + 1) - 1, n + 1)
  counts <- vapply(trees, function(tree) {
    ed <- tree$edge
    min(CladeSizes(tree, nodes = ed[ed[, 1] == n + 3, 2]))
  }, double(1)) |>
    table() |>
    as.numeric()
  expect_equal(exp(c(.LogR(1, n))) * NRooted(n), counts)
  
})

test_that(".LogD() succeeds", {
  # D is the probability that, in a randomly selected tree on `leaves`, the
  # smaller subclade of m taxa will receive taxa with labels `drawn`
  expect_equal(.LogD(c(4, 2), c(4, 2)), 0)
  expect_equal(.LogD(c(4, 0), c(4, 2)), 0 - lchoose(6, 4))
  expect_equal(.LogD(c(2, 2), c(4, 2)), lchoose(4, 2) - lchoose(6, 4))
  expect_equal(.LogD(c(2, 2, 2), c(4, 2, 2)), lchoose(4, 2) - lchoose(8, 6))
  expect_equal(.LogD(c(2, 2, 1), c(4, 2, 2)),
               lchoose(4, 2) + lchoose(2, 1) - lchoose(8, 5))
})

test_that(".LogD() succeeds - balanced trees", {
  # As defined, the 'smaller' subclade is drawn at random if both are the same
  # size
  expect_equal(.LogD(c(1, 0), c(1, 1)), log(1 / 2))
  expect_equal(.LogD(c(2, 0), c(2, 2)), log(1/3 / 2))
})

test_that(".ValidDraws() succeeds", {
  # Drawing the contents of the smallest clade.
  expect_equal(unname(.ValidDraws(c(0, 4, 0))), cbind(c(0, 0), 1:2, c(0, 0)))
  
  # As 2, 0, 0 == 0, 2, 0, only one of these combinations should be listed
  expect_equal(unname(.ValidDraws(c(2, 2, 0))),
               cbind(c(1, 0, 1, 0), c(0, 1, 1, 2), rep(0, 4)))
})

test_that(".LogB() solves simple cases", {
  # B(b | tokens) is the probability that state b is reconstructed at the base
  # of a clade with leaves labelled `tokens`
  expect_error(.LogB(0, c(0, 1, 0)), "token. must be 1..")
  
  # Singletons - special cases
  expect_equal(.LogB(1, c(0, 1, 0)), log(0))
  expect_equal(.LogB(2, c(0, 1, 0)), log(1))
  expect_equal(.LogB(3, c(0, 1, 0)), log(0))
  
  # Pairs - special cases
  expect_equal(.LogB(1, c(0, 2, 0)), log(0))
  expect_equal(.LogB(2, c(0, 2, 0)), log(1))
  expect_equal(.LogB(3, c(0, 2, 0)), log(0))
  
  expect_equal(.LogB(1, c(1, 1, 0)), log(0))
  expect_equal(.LogB(2, c(1, 1, 0)), log(0))
  expect_equal(.LogB(3, c(1, 1, 0)), log(1))
  
  # Simpletons
  expect_equal(.LogB(1, c(0, 10, 0)), log(0))
  expect_equal(.LogB(2, c(0, 10, 0)), log(1))
  expect_equal(.LogB(3, c(0, 10, 0)), log(0))
})

test_that(".LogB() arithmetic is correct", {
  expect_equal(exp(LogSumExp(.LogB(1, c(2, 1, 0)),
                             .LogB(2, c(2, 1, 0)),
                             .LogB(3, c(2, 1, 0)))), 1)
  
  expect_equal(exp(.LogB(1, c(1, 2, 0))) * 3, 0)
  expect_equal(exp(.LogB(2, c(1, 2, 0))) * 3, 2) # (2a, (2b, 1)), (2b, (2a, 1))
  expect_equal(exp(.LogB(3, c(1, 2, 0))) * 3, 1) # (1, (2, 2))
  
  # 15 rooted four-leaf trees
  trees <- as.phylo(0:14, 5)
  if (interactive()) {
    dev.new()
    par(mfrow = c(5, 3), cex = 0.9, mar = rep(0.4, 4), xpd = NA)
    on.exit(dev.off())
  }
  rootState <- (vapply(trees, PlotCharacter, matrix(NA, 9, 2),
                       StringToPhyDat("?1122"), 1, plot = interactive())[6, , ] *
                  1:2) |>
    colSums() |> 
    table()

  expect_equal(exp(.LogB(1, c(2, 2, 0))) * 15, rootState[[1]]) # 2
  expect_equal(exp(.LogB(2, c(2, 2, 0))) * 15, rootState[[2]]) # 2
  expect_equal(exp(.LogB(3, c(2, 2, 0))) * 15, rootState[[3]]) # 11
  
  
  expect_equal(exp(LogSumExp(.LogB(1, c(6, 3, 0)),
                             .LogB(2, c(6, 3, 0)),
                             .LogB(3, c(6, 3, 0)))), 1)
  expect_equal(exp(LogSumExp(.LogB(1, c(4, 3, 2)),
                             .LogB(2, c(4, 3, 2)),
                             .LogB(3, c(4, 3, 2)))), 1)
})



test_that("MaddisonSlatkin() recursion bottoms", {
  expect_equal(MaddisonSlatkin(0, c(1, 1)), log(0))
  
  expect_equal(MaddisonSlatkin(1, c(1, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(2, 0)), log(1))
  expect_equal(MaddisonSlatkin(1, c(1, 0, 0, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(0, 0, 0, 2)), log(1))
})
  
  
test_that("MaddisonSlatkin() is numerically correct", {
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
  
  exp(MaddisonSlatkin(2, c(2, 2)))
  exp(MaddisonSlatkin(1, c(2, 2)))
  LogCarter1(1,2,2)
  LogCarter1(2,2,2)
  
  
  expect_equal(MaddisonSlatkin(1, c(8, 24)) - LnRooted(32),
               LogCarter1(1, 8, 24))
  LogCarter1(2, 8, 24)
  LogCarter1(3, 8, 24)
})
