test_that("Steel 1996 formulae are useful", {
  library("TreeTools", quietly = TRUE)
  
  nTip <- 7
  tree <- BalancedTree(nTip)
  dist <- c(2, 3, 2) / nTip
  kMax <- 4
  
  logProbs <- active_parsimony_dist(tree, dist, kMax)
  calcProbs <- exp(logProbs)
  expect_equal(sum(calcProbs), 1)
  
  set.seed(1)
  nSim <- 2000
  simTable <- vapply(
    replicate(nSim, sample(seq_along(dist), nTip, replace = TRUE, prob = dist) |>
                cbind() |>
                `rownames<-`(TipLabels(nTip)) |>
                MatrixToPhyDat(), simplify = FALSE),
    function(char) TreeLength(tree, char), 1) |>
    table()
  
  simCount <- setNames(rep(0, length(calcProbs)), names(calcProbs))
  commonNames <- intersect(names(simTable), names(calcProbs))
  simCount[commonNames] <- as.numeric(simTable[commonNames])
  
  chi2 <- suppressWarnings(chisq.test(x = simCount, p = calcProbs))
  expect_gt(chi2$p.value, 0.001) # i.e. NSD
  
  expMean <- sum(as.numeric(names(calcProbs)) * calcProbs)
  simMean <- sum(as.numeric(names(simTable)) * simTable) / sum(simTable)
  expect_equal(simMean, expMean, tolerance = 0.1)
})
