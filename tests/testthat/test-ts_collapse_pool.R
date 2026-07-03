# ts_collapse_pool is an internal (unexported) C++ bridge reachable via
# TreeSearch:::ts_collapse_pool.  It assumes every edge matrix it is fed is
# fully binary (the engine's own pool is always binary); a caller who
# hand-builds a non-binary edge matrix and bypasses MaximizeParsimony()'s own
# RootTree/MakeTreeBinary normalization can otherwise hang the C++ reroot/
# postorder walk forever (T-332).

test_that("ts_collapse_pool errors cleanly (and quickly) on a non-binary tree", {
  charDat <- StringToPhyDat(c("0000", "0011", "0101", "1010"), c("a", "b", "c", "d"))
  ds <- make_ts_data(charDat)
  scoringConfig <- list(
    min_steps = integer(0), concavity = Inf, xpiwe = FALSE,
    xpiwe_r = 0.5, xpiwe_max_f = 5.0, obs_count = integer(0),
    infoAmounts = NULL
  )

  # Star tree: node 5 has four children (a, b, c, d) -- not binary.
  starEdge <- cbind(c(5L, 5L, 5L, 5L), c(1L, 2L, 3L, 4L))

  elapsed <- system.time(
    result <- tryCatch(
      TreeSearch:::ts_collapse_pool(
        list(starEdge), ds$contrast, ds$tip_data, ds$weight, ds$levels,
        scoringConfig
      ),
      error = function(e) e
    )
  )["elapsed"]

  expect_true(inherits(result, "error"))
  expect_match(conditionMessage(result), "not binary", ignore.case = TRUE)
  # Prior to the T-332 fix this call hung indefinitely (killed at 15s); a
  # generous ceiling here proves termination without being a timing test.
  expect_lt(elapsed, 10)
})

test_that("ts_collapse_pool succeeds on a genuinely binary tree", {
  charDat <- StringToPhyDat(c("0000", "0011", "0101", "1010"), c("a", "b", "c", "d"))
  ds <- make_ts_data(charDat)
  scoringConfig <- list(
    min_steps = integer(0), concavity = Inf, xpiwe = FALSE,
    xpiwe_r = 0.5, xpiwe_max_f = 5.0, obs_count = integer(0),
    infoAmounts = NULL
  )

  # Binary rooted tree: ((a,b),(c,d)); root = 5, internals = 6, 7.
  binaryEdge <- cbind(
    c(5L, 5L, 6L, 6L, 7L, 7L),
    c(6L, 7L, 1L, 2L, 3L, 4L)
  )

  result <- TreeSearch:::ts_collapse_pool(
    list(binaryEdge), ds$contrast, ds$tip_data, ds$weight, ds$levels,
    scoringConfig
  )

  expect_true(is.list(result))
  expect_true(result$n_topologies >= 1L)
})
