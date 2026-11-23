# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Test cross-replicate consensus constraint tightening.
# When consensusConstrain = TRUE, the strict consensus of pool trees is
# enforced as topological constraints for subsequent replicates, focusing
# search on uncertain parts of the tree.

ts_driven <- function(ds, maxReps = 10L, targetHits = 3L,
                      rssRounds = 1L, ratchetCycles = 3L,
                      driftCycles = 0L, xssRounds = 1L,
                      xssPartitions = 4L, verbosity = 0L,
                      consensusConstrain = FALSE) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReps, targetHits = targetHits,
    tbrMaxHits = 1L, ratchetCycles = ratchetCycles,
    ratchetPerturbProb = 0.04, ratchetPerturbMode = 0L,
    ratchetPerturbMaxMoves = 0L, ratchetAdaptive = FALSE,
    driftCycles = driftCycles, driftAfdLimit = 3L,
    driftRfdLimit = 0.1,
    xssRounds = xssRounds, xssPartitions = xssPartitions,
    rssRounds = rssRounds, cssRounds = 0L, cssPartitions = 4L,
    sectorMinSize = 6L, sectorMaxSize = 50L,
    fuseInterval = 3L, fuseAcceptEqual = FALSE,
    poolMaxSize = 100L, poolSuboptimal = 0,
    maxSeconds = 0, verbosity = verbosity,
    nThreads = 1L, tabuSize = 100L,
    sprFirst = FALSE, wagnerStarts = 1L,
    consensusStableReps = 0L, adaptiveLevel = FALSE,
    consensusConstrain = consensusConstrain
  )
}

data("inapplicable.phyData", package = "TreeSearch")
phy <- inapplicable.phyData[["Vinther2008"]]
ds <- make_ts_data(phy)

test_that("consensus constraint tightening produces valid trees", {
  set.seed(7129)
  result <- ts_driven(ds, maxReps = 10L, targetHits = 6L,
                      ratchetCycles = 3L, xssRounds = 1L,
                      consensusConstrain = TRUE)
  expect_true(result$best_score > 0)
  expect_true(result$replicates >= 1L)
  expect_true(result$pool_size >= 1L)
  n_tip <- length(phy)
  edges <- result$trees[[1]]
  expect_equal(nrow(edges), 2L * (n_tip - 1L))
})

test_that("consensus constraint doesn't degrade score quality", {
  set.seed(4163)
  without <- ts_driven(ds, maxReps = 10L, targetHits = 5L,
                       ratchetCycles = 3L, xssRounds = 1L,
                       consensusConstrain = FALSE)

  set.seed(4163)
  with_cc <- ts_driven(ds, maxReps = 10L, targetHits = 5L,
                       ratchetCycles = 3L, xssRounds = 1L,
                       consensusConstrain = TRUE)

  # Constraint tightening should not produce a worse score
  expect_true(with_cc$best_score <= without$best_score + 1)
})

# Test on larger dataset where constraint tightening can make a difference
phy_big <- inapplicable.phyData[["Agnarsson2004"]]
ds_big <- make_ts_data(phy_big)

test_that("consensus constraint works on larger dataset", {
  set.seed(6842)
  result <- ts_driven(ds_big, maxReps = 8L, targetHits = 4L,
                      ratchetCycles = 3L, driftCycles = 1L,
                      xssRounds = 2L,
                      consensusConstrain = TRUE)
  expect_true(result$best_score > 0)
  expect_true(result$pool_size >= 1L)
  n_tip <- length(phy_big)
  edges <- result$trees[[1]]
  expect_equal(nrow(edges), 2L * (n_tip - 1L))
})
