# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for consensus-stability stopping and adaptive search level.
# These features require multiple replicates, so tests are intentionally
# lightweight (small datasets, low replicate caps).

ts_driven <- TreeSearch:::ts_driven_search

test_that("consensusStableReps = 0 does not trigger early stop", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])

  set.seed(7832)
  result <- ts_driven(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 6L, targetHits = 20L,  # high target so hits won't fire
    ratchetCycles = 2L, driftCycles = 0L,
    xssRounds = 0L, rssRounds = 0L, cssRounds = 0L,
    consensusStableReps = 0L, adaptiveLevel = FALSE,
    verbosity = 0L
  )
  # Should run all 6 replicates (consensus disabled, target unreachable)
  expect_equal(result$replicates, 6L)
  expect_false(result$consensus_stable)
})

test_that("consensus stability stops search early", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])

  set.seed(2491)
  result <- ts_driven(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 30L, targetHits = 25L,  # high target to let consensus fire
    ratchetCycles = 3L, driftCycles = 1L,
    xssRounds = 1L, rssRounds = 0L, cssRounds = 0L,
    consensusStableReps = 3L, adaptiveLevel = FALSE,
    verbosity = 0L
  )
  # Should stop before max replicates
  expect_lt(result$replicates, 30L)
  expect_true(result$consensus_stable)
  expect_gt(length(result$trees), 0L)
})

test_that("adaptive level adjusts without crashing", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])

  set.seed(6153)
  result <- ts_driven(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 8L, targetHits = 6L,
    ratchetCycles = 3L, driftCycles = 1L,
    xssRounds = 1L, rssRounds = 0L, cssRounds = 0L,
    consensusStableReps = 0L, adaptiveLevel = TRUE,
    verbosity = 0L
  )
  expect_gt(result$replicates, 0L)
  expect_gt(length(result$trees), 0L)
  expect_false(result$consensus_stable)
})

test_that("both features work together", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])

  set.seed(3847)
  result <- ts_driven(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 30L, targetHits = 25L,
    ratchetCycles = 3L, driftCycles = 1L,
    xssRounds = 1L, rssRounds = 0L, cssRounds = 0L,
    consensusStableReps = 3L, adaptiveLevel = TRUE,
    verbosity = 0L
  )
  expect_gt(result$replicates, 0L)
  expect_gt(length(result$trees), 0L)
  # Either consensus stability or max replicates
  expect_true(result$consensus_stable || result$replicates == 30L)
})

test_that("SearchControl round-trips new parameters", {
  sc <- SearchControl(consensusStableReps = 5L, adaptiveLevel = TRUE)
  expect_equal(sc$consensusStableReps, 5L)
  expect_true(sc$adaptiveLevel)
  expect_s3_class(sc, "SearchControl")

  # Defaults
  sc0 <- SearchControl()
  expect_equal(sc0$consensusStableReps, 0L)
  expect_false(sc0$adaptiveLevel)
})
