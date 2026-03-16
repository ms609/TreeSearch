library("TreeTools")

# ===========================================================================
# Tests for NA-aware incremental scoring (Phase 2A)
# ===========================================================================

make_ts_data <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  weight <- at$weight
  levels <- at$levels
  list(contrast = contrast, tip_data = tip_data,
       weight = weight, levels = levels)
}

# Helper: run driven search with sensible defaults for testing
ts_driven <- function(ds, maxReplicates = 3L, targetHits = 2L,
                      ratchetCycles = 3L, xssRounds = 1L,
                      xssPartitions = 2L, fuseInterval = 2L,
                      driftCycles = 2L,
                      maxSeconds = 0, verbosity = 0L,
                      concavity = -1.0, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    ratchetCycles = ratchetCycles,
    xssRounds = xssRounds,
    xssPartitions = xssPartitions,
    fuseInterval = fuseInterval,
    driftCycles = driftCycles,
    maxSeconds = maxSeconds,
    verbosity = verbosity,
    concavity = concavity,
    ...
  )
}

# --- Test 1: Driven search on inapplicable datasets ---

test_that("Driven search on inapplicable datasets finds good scores", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)

  # Upper bounds are generous — these are quick searches with few replicates.
  # The goal is to verify the search runs and improves, not find optimal.
  test_cases <- list(
    list("Vinther2008", 120),
    list("Agnarsson2004", 900),
    list("Aguado2009", 650)
  )

  for (tc in test_cases) {
    ds_name <- tc[[1]]
    upper_bound <- tc[[2]]
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)

    set.seed(7342)
    result <- ts_driven(ds, maxReplicates = 5L, targetHits = 3L,
                        maxSeconds = 30)

    expect_true(result$best_score <= upper_bound,
      label = paste(ds_name, "score", result$best_score, "<=", upper_bound))
    expect_true(result$best_score > 0,
      label = paste(ds_name, "score is positive"))
  }
})


# --- Test 2: set.seed reproducibility on inapplicable datasets ---

test_that("Driven search on NA data is reproducible with set.seed", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  run_search <- function() {
    set.seed(8819)
    ts_driven(ds, maxReplicates = 3L, targetHits = 2L)
  }

  r1 <- run_search()
  r2 <- run_search()
  expect_equal(r1$best_score, r2$best_score)
  expect_equal(r1$trees, r2$trees)
})


# --- Test 3: NA datasets with implied weights ---

test_that("Driven search on NA data with implied weights", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  set.seed(3156)
  result <- ts_driven(ds, maxReplicates = 3L, targetHits = 2L,
                      concavity = 10.0)

  expect_true(result$best_score > 0)
  expect_true(result$best_score < Inf)
  expect_true(length(result$trees) >= 1)
})


# --- Test 4: Multiple inapplicable datasets ---

test_that("Driven search on multiple inapplicable datasets", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)

  ds_names <- c("Vinther2008", "Aguado2009", "Agnarsson2004",
                "Aria2015", "DeAssis2011")

  for (ds_name in ds_names) {
    if (is.null(inapplicable.phyData[[ds_name]])) next
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)

    set.seed(5523)
    result <- ts_driven(ds, maxReplicates = 2L, targetHits = 2L,
                        maxSeconds = 20)

    expect_true(result$best_score > 0,
      label = paste(ds_name, ": positive score"))
    expect_true(result$best_score < Inf,
      label = paste(ds_name, ": finite score"))
    expect_true(length(result$trees) >= 1,
      label = paste(ds_name, ": at least one tree"))
  }
})


# --- Test 5: Score verification with ts_fitch_score ---

test_that("Returned trees from NA search have correct scores", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  set.seed(6091)
  result <- ts_driven(ds, maxReplicates = 3L, targetHits = 2L)

  for (i in seq_along(result$trees)) {
    edge_i <- result$trees[[i]]
    verify_score <- TreeSearch:::ts_fitch_score(
      edge_i, ds$contrast, ds$tip_data, ds$weight, ds$levels,
      concavity = -1.0)
    expect_equal(verify_score, result$scores[i], tolerance = 1e-10,
      label = paste("Tree", i, "score matches"))
  }
})


# --- Test 6: Pool collects suboptimal trees on NA datasets ---

test_that("Pool collects suboptimal trees on NA data", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  set.seed(2487)
  result <- ts_driven(ds, maxReplicates = 5L, targetHits = 3L,
                      poolSuboptimal = 5.0, poolMaxSize = 20L)

  expect_true(length(result$trees) >= 1)
  expect_true(all(result$scores <= result$best_score + 5.0 + 1e-10))
})


# --- Test 7: Timeout works on NA datasets ---

test_that("Timeout works on NA datasets", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  set.seed(1103)
  result <- ts_driven(ds, maxReplicates = 1000L, targetHits = 1000L,
                      maxSeconds = 0.5)

  expect_true(result$timed_out)
  expect_true(result$best_score > 0)
})


# --- Test 8: NA + IW score verification ---

test_that("NA + IW returned scores are verified correct", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  k <- 10.0

  set.seed(9471)
  result <- ts_driven(ds, maxReplicates = 3L, targetHits = 2L,
                      concavity = k)

  for (i in seq_along(result$trees)) {
    verify <- TreeSearch:::ts_fitch_score(
      result$trees[[i]], ds$contrast, ds$tip_data, ds$weight, ds$levels,
      concavity = k)
    expect_equal(verify, result$scores[i], tolerance = 1e-8,
      label = paste("IW tree", i, "score matches"))
  }
})


# --- Test 9: Score verification on NA datasets with more replicates ---

test_that("More replicates improve NA search scores", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  # Short search
  set.seed(4419)
  short <- ts_driven(ds, maxReplicates = 1L, targetHits = 1L,
                     ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L)

  # Longer search
  set.seed(4419)
  long <- ts_driven(ds, maxReplicates = 5L, targetHits = 3L)

  # Longer should be at least as good
  expect_true(long$best_score <= short$best_score)
})
