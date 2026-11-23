# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Test that conflict-guided sector selection in RSS works correctly.
# The feature biases random sector selection toward tree regions that are
# uncertain across pool trees (splits not present in all best-score trees).
#
# We test end-to-end through driven_search, which computes the split
# frequency table from the pool and passes it to run_single_replicate → RSS.

ts_driven <- function(ds, maxReps = 10L, targetHits = 3L,
                      rssRounds = 1L, ratchetCycles = 3L,
                      driftCycles = 0L, xssRounds = 1L,
                      xssPartitions = 4L, nThreads = 1L,
                      verbosity = 0L) {
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
    nThreads = nThreads, tabuSize = 100L,
    sprFirst = FALSE, wagnerStarts = 1L,
    consensusStableReps = 0L, adaptiveLevel = FALSE
  )
}

# Use Vinther2008 (23 tips) — small enough for fast testing,
# large enough for RSS sectors.
data("inapplicable.phyData", package = "TreeSearch")
phy <- inapplicable.phyData[["Vinther2008"]]
ds <- make_ts_data(phy)

test_that("driven search with RSS produces valid results", {
  set.seed(8412)
  # Multiple replicates so the pool accumulates >1 tree,

  # triggering conflict-guided selection in later replicates.
  result <- ts_driven(ds, maxReps = 8L, targetHits = 5L,
                      rssRounds = 1L, ratchetCycles = 3L,
                      xssRounds = 1L)
  expect_true(result$best_score > 0)
  expect_true(result$replicates >= 1L)
  expect_true(result$pool_size >= 1L)
  # Validate at least one tree in the pool
  edges <- result$trees[[1]]
  n_tip <- length(phy)
  expect_equal(nrow(edges), 2L * (n_tip - 1L))
})

test_that("RSS conflict path doesn't degrade score quality", {
  set.seed(2917)
  # Run without RSS
  no_rss <- ts_driven(ds, maxReps = 6L, targetHits = 4L,
                      rssRounds = 0L, ratchetCycles = 3L,
                      xssRounds = 1L)

  set.seed(2917)
  # Same seed, with RSS (conflict guidance active after pool fills)
  with_rss <- ts_driven(ds, maxReps = 6L, targetHits = 4L,
                        rssRounds = 1L, ratchetCycles = 3L,
                        xssRounds = 1L)

  # Score with RSS should be no worse (allowing for stochastic variation)
  # We don't require strict improvement — just no regression.
  expect_true(with_rss$best_score <= no_rss$best_score + 1)
})

# Use a larger dataset where RSS has more room to work.
phy_big <- inapplicable.phyData[["Agnarsson2004"]]
ds_big <- make_ts_data(phy_big)

test_that("conflict-guided RSS works on larger dataset", {
  set.seed(5039)
  result <- ts_driven(ds_big, maxReps = 5L, targetHits = 3L,
                      rssRounds = 2L, ratchetCycles = 3L,
                      driftCycles = 1L, xssRounds = 2L)
  expect_true(result$best_score > 0)
  expect_true(result$pool_size >= 1L)
  n_tip <- length(phy_big)
  edges <- result$trees[[1]]
  expect_equal(nrow(edges), 2L * (n_tip - 1L))
})
