# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for parallel driven search (Phase 5).
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

data("inapplicable.phyData", package = "TreeSearch")
vinther <- inapplicable.phyData[["Vinther2008"]]

# --- 1. Serial equivalence (nThreads=1) ---

test_that("nThreads=1 produces identical results to default serial search", {
  skip_on_cran()
  set.seed(8741)
  result_default <- TreeSearch:::ts_driven_search(
    contrast = attributes(vinther)$contrast,
    tip_data = matrix(unlist(vinther, use.names = FALSE),
                      nrow = length(vinther), byrow = TRUE),
    weight = attributes(vinther)$weight,
    levels = attributes(vinther)$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    nThreads = 1L
  )

  set.seed(8741)
  result_serial <- TreeSearch:::ts_driven_search(
    contrast = attributes(vinther)$contrast,
    tip_data = matrix(unlist(vinther, use.names = FALSE),
                      nrow = length(vinther), byrow = TRUE),
    weight = attributes(vinther)$weight,
    levels = attributes(vinther)$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )

  expect_equal(result_default$best_score, result_serial$best_score)
  expect_equal(result_default$replicates, result_serial$replicates)
})

# --- 2. Parallel correctness ---

test_that("Parallel search (2 threads) produces valid trees with correct scores", {
  skip_on_cran()
  set.seed(3192)
  ds <- make_ts_data(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)
  expect_true(result$pool_size > 0)
  expect_true(is.finite(result$best_score))

  # Verify each tree is valid and scores correctly
  for (i in seq_along(result$trees)) {
    edge <- result$trees[[i]]
    score <- TreeSearch:::ts_fitch_score(
      edge, ds$contrast, ds$tip_data, ds$weight, ds$levels
    )
    expect_equal(score, result$scores[i], tolerance = 0.01,
                 info = paste("Tree", i, "score mismatch"))
    n_tip <- length(vinther)
    expect_equal(nrow(edge), 2 * (n_tip - 1))
  }
})

# --- 3. Timeout in parallel mode ---

test_that("Parallel search respects timeout", {
  skip_on_cran()
  # Use a larger dataset so each replicate takes long enough to trigger
  # timeout reliably on fast hardware (23-tip Vinther completes <1ms/rep)
  agnarsson <- inapplicable.phyData[["Agnarsson2004"]]
  ds_lg <- make_ts_data(agnarsson)
  t0 <- proc.time()["elapsed"]
  result <- TreeSearch:::ts_driven_search(
    contrast = ds_lg$contrast, tip_data = ds_lg$tip_data,
    weight = ds_lg$weight, levels = ds_lg$levels,
    maxReplicates = 1000L, targetHits = 999L,
    maxSeconds = 2.0, verbosity = 0L,
    nThreads = 2L
  )
  elapsed <- proc.time()["elapsed"] - t0

  # `timed_out` is set only on the deadline path (never on natural completion),
  # so this alone proves the timeout fired.
  expect_true(result$timed_out)
  # And the search actually stopped rather than running the full budget.  This
  # is a working-vs-broken discriminator, NOT a wall-clock precision assertion:
  # a working timeout returns in seconds-to-low-tens (bounded by one in-flight
  # sectorial pass per worker), whereas a broken one runs all 1000 replicates
  # (>1000 s) — a ~100x gap, so 60 s absorbs CI/shared-VM jitter while still
  # catching an unbounded run.  The old 15 s bound was too tight: the per-test
  # wall time here is observed to vary from ~5 s to ~30 s on the same machine.
  expect_lt(elapsed, 60)
  expect_lt(result$replicates, 1000L)  # stopped before exhausting the budget
})

# --- 4. Edge cases ---

test_that("nThreads=0 (auto-detect) is capped at 2 in tests", {
  skip_on_cran()
  # nThreads=0 auto-detects CPU count; use nThreads=2 to stay within

  # the 2-core-per-agent limit (see AGENTS.md).
  ds <- make_ts_data(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)
})

test_that("nThreads > maxReplicates is clamped", {
  skip_on_cran()
  ds <- make_ts_data(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L,
    nThreads = 100L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$replicates <= 2L)
})

test_that("Single replicate in parallel mode works", {
  skip_on_cran()
  ds <- make_ts_data(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 1L, targetHits = 1L, verbosity = 0L,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_equal(result$replicates, 1L)
})

# --- 5. IW + parallel ---

test_that("Implied weights works in parallel mode", {
  skip_on_cran()
  ds <- make_ts_data(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    concavity = 10.0,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)

  # Verify IW scores
  for (i in seq_along(result$trees)) {
    score <- TreeSearch:::ts_fitch_score(
      result$trees[[i]], ds$contrast, ds$tip_data,
      ds$weight, ds$levels, concavity = 10.0
    )
    expect_equal(score, result$scores[i], tolerance = 0.01)
  }
})

# --- 6. NA + parallel ---

test_that("Inapplicable characters work in parallel mode", {
  skip_on_cran()
  ds <- make_ts_data(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    nThreads = 2L
  )

  # Vinther2008 has inapplicable characters
  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)
})

# --- 7. R-level MaximizeParsimony with parallel ---

test_that("MaximizeParsimony with nThreads > 1 works end-to-end", {
  skip_on_cran()
  set.seed(5023)
  result <- MaximizeParsimony(vinther, maxReplicates = 3L,
                              targetHits = 2L, nThreads = 2L,
                              verbosity = 0L)

  expect_s3_class(result, "multiPhylo")
  expect_true(length(result) > 0)
  expect_true(attr(result, "score") > 0)
  expect_true(is.finite(attr(result, "score")))

  # All trees should have correct number of tips
  for (i in seq_along(result)) {
    expect_equal(length(result[[i]]$tip.label), length(vinther))
    expect_true(ape::is.binary(result[[i]]))
  }
})

# --- 8. Parallel hits_to_best matches serial ---

test_that("Parallel hits_to_best tracks independent replicate hits", {
  skip_on_cran()
  # Bug (T-242): extract_into() rebuilt hits_to_best from pool entries,

  # losing the real count.  A 1-topology pool always reported 1 hit.
  agn <- inapplicable.phyData[["Agnarsson2004"]]
  ds <- make_ts_data(agn)

  set.seed(6291)
  r_serial <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 15L, targetHits = 15L,
    verbosity = 0L, nThreads = 1L
  )

  set.seed(6291)
  r_par <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 15L, targetHits = 15L,
    verbosity = 0L, nThreads = 2L
  )

  # Parallel hits should be in the same ballpark as serial, not 1
  expect_true(r_par$hits_to_best >= 2,
              info = paste("Parallel hits:", r_par$hits_to_best,
                           "Serial hits:", r_serial$hits_to_best))
  # And serial should also have multiple hits on this dataset
  expect_true(r_serial$hits_to_best >= 2)
})

# --- 9. Pool suboptimal in parallel ---


test_that("Pool suboptimal collection works in parallel", {
  skip_on_cran()
  ds <- make_ts_data(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    poolSuboptimal = 5.0,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  # With suboptimal > 0, we might get multiple trees
  # All scores should be within suboptimal of best
  for (i in seq_along(result$scores)) {
    expect_true(result$scores[i] <= result$best_score + 5.0 + 0.01)
  }
})
