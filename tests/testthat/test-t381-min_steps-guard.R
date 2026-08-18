# T-381: min_steps invariant violation guard
# Test that an invalid min_steps (where precomputed_steps > min_steps_r)
# raises an error rather than silently clamping to 0.

skip_on_cran()

test_that("T-381: min_steps < precomputed_steps raises error", {
  # Create a test dataset with autapomorphies that will have precomputed_steps > 0.
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,  # informative binary
    2, 0, 0, 0, 0, 0, 0, 0   # autapomorphy in tip 1 -> will have precomputed_steps = 1
  ), nrow = 8, byrow = TRUE, dimnames = list(paste0("t", 1:8), NULL))

  dataset <- TreeTools::MatrixToPhyDat(mat)
  attrs <- attributes(dataset)

  # Broken min_steps: set min_steps[1] = 0, but char 2 will have precomputed_steps = 1
  # This violates the invariant: min_steps_r[1] (0) < precomputed_steps[1] (1)
  broken_min_steps <- c(10L, 0L)

  # Call ts_resample_search with invalid min_steps.
  # This triggers build_dataset internally and should raise an error
  # instead of silently mis-scoring.
  expect_error(
    TreeSearch:::ts_resample_search(
      contrast = attrs$contrast,
      tip_data = matrix(as.integer(unlist(dataset, use.names = FALSE)),
                        nrow = length(dataset)),
      weight = attrs$weight,
      levels = attrs$levels,
      maxReplicates = 1L, targetHits = 1L,
      min_steps = broken_min_steps
    ),
    "Internal invariant violation.*precomputed_steps.*min_steps"
  )
})

test_that("T-381: valid min_steps passes", {
  # Sanity check: with proper min_steps (>= precomputed_steps),
  # ts_resample_search should work.
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1
  ), nrow = 8, byrow = TRUE, dimnames = list(paste0("t", 1:8), NULL))

  dataset <- TreeTools::MatrixToPhyDat(mat)
  attrs <- attributes(dataset)

  # Valid min_steps: large enough that precomputed_steps <= min_steps
  valid_min_steps <- c(100L)

  result <- TreeSearch:::ts_resample_search(
    contrast = attrs$contrast,
    tip_data = matrix(as.integer(unlist(dataset, use.names = FALSE)),
                      nrow = length(dataset)),
    weight = attrs$weight,
    levels = attrs$levels,
    maxReplicates = 1L, targetHits = 1L,
    min_steps = valid_min_steps
  )

  # Should complete without error and return a valid result
  expect_true(is.numeric(result$score))
  expect_gte(result$score, 0)
})
