# Extracted from test-ts-simd.R:144

# prequel ----------------------------------------------------------------------
morphy_ew_ref <- function(tree, dataset) {
  suppressWarnings(TreeSearch::Fitch(tree, dataset))
}

# test -------------------------------------------------------------------------
dataset <- inapplicable.phyData[["Vinther2008"]]
ds <- make_ts_data(dataset)
set.seed(8371)
result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data,
    ds$weight, ds$levels,
    maxReplicates = 2L, targetHits = 1L,
    ratchetCycles = 1L, driftCycles = 0L,
    xssPartitions = 2L, rssRounds = 0L, cssRounds = 0L,
    cssPartitions = 2L, fuseInterval = 0L,
    poolMaxSize = 2L, poolSuboptimal = 0,
    ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 50L,
    ratchetAdaptive = FALSE, maxSeconds = 30,
    verbosity = 0L
  )
