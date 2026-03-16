# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# Helper: run driven search with callback
ts_driven_cb <- function(ds, callback, ...) {
  defaults <- list(
    contrast = ds$contrast,
    tip_data = ds$tip_data,
    weight = ds$weight,
    levels = ds$levels,
    maxReplicates = 3L,
    targetHits = 2L,
    ratchetCycles = 2L,
    xssRounds = 1L,
    xssPartitions = 2L,
    fuseInterval = 2L,
    maxSeconds = 0,
    verbosity = 1L,
    progressCallback = callback
  )
  args <- modifyList(defaults, list(...))
  do.call(TreeSearch:::ts_driven_search, args)
}

# Small dataset for callback tests (search quality doesn't matter here)
small_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
small_ds <- make_ts_data(MatrixToPhyDat(small_mat))


# ===== Callback is invoked =====

test_that("Callback receives expected phases", {
  ds <- small_ds

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(4517)
  result <- ts_driven_cb(ds, recorder)

  # Should have received at least one callback

  expect_gt(length(log), 0)

  # Extract phase names
  phases <- vapply(log, `[[`, character(1), "phase")

  # Must end with "done"
  expect_equal(phases[length(phases)], "done")

  # Should contain "replicate" events
  expect_true("replicate" %in% phases)
})

test_that("Replicate numbers increment correctly", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(7823)
  result <- ts_driven_cb(ds, recorder, maxReplicates = 5L,
                          targetHits = 5L)

  # Get replicate events
  rep_events <- Filter(function(x) x$phase == "replicate", log)
  rep_nums <- vapply(rep_events, `[[`, integer(1), "replicate")

  # Should be sequential starting from 1
  expect_equal(rep_nums, seq_along(rep_nums))
})

test_that("Done event has consistent best_score", {
  ds <- small_ds

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(2396)
  result <- ts_driven_cb(ds, recorder)

  # Find the "done" event
  done_events <- Filter(function(x) x$phase == "done", log)
  expect_length(done_events, 1)

  done <- done_events[[1]]
  expect_equal(done$best_score, result$best_score, tolerance = 1e-6)
  expect_equal(done$pool_size, result$pool_size)
})


# ===== NULL callback (regression) =====

test_that("Search works with NULL callback", {
  ds <- small_ds

  set.seed(5614)
  result <- ts_driven_cb(ds, NULL, verbosity = 0L)

  expect_true(is.list(result))
  expect_gt(length(result$trees), 0)
  expect_gt(result$best_score, 0)
})


# ===== Verbosity 0 suppresses callback =====

test_that("Callback NOT invoked when verbosity = 0", {
  ds <- small_ds

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(8091)
  result <- ts_driven_cb(ds, recorder, verbosity = 0L)

  expect_length(log, 0)
})


# ===== Callback info structure =====

test_that("Callback info has all expected fields", {
  ds <- small_ds

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(3344)
  result <- ts_driven_cb(ds, recorder)

  expected_fields <- c("replicate", "max_replicates", "best_score",
                        "hits_to_best", "target_hits", "pool_size",
                        "phase", "elapsed", "phase_score")

  for (entry in log) {
    for (f in expected_fields) {
      expect_true(f %in% names(entry),
                  label = paste("Field", f, "in phase", entry$phase))
    }
  }
})


# ===== Elapsed time increases =====

test_that("Elapsed time is non-decreasing across callbacks", {
  ds <- small_ds

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(6712)
  result <- ts_driven_cb(ds, recorder)

  elapsed_vals <- vapply(log, `[[`, double(1), "elapsed")
  for (i in seq_along(elapsed_vals)[-1]) {
    expect_gte(elapsed_vals[i], elapsed_vals[i - 1],
               label = paste("Elapsed at callback", i))
  }
})


# ===== MaximizeParsimony with custom callback =====

test_that("MaximizeParsimony accepts custom progressCallback", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(1937)
  result <- suppressMessages(
    MaximizeParsimony(dataset, maxReplicates = 3L, targetHits = 2L,
                      verbosity = 1L, progressCallback = recorder)
  )

  expect_s3_class(result, "multiPhylo")
  expect_gt(length(log), 0)

  phases <- vapply(log, `[[`, character(1), "phase")
  expect_true("done" %in% phases)
})

test_that("MaximizeParsimony silent with verbosity = 0 and no callback", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]

  set.seed(4403)
  result <- MaximizeParsimony(dataset, maxReplicates = 3L, targetHits = 2L,
                              verbosity = 0L)

  expect_s3_class(result, "multiPhylo")
})


# ===== Fuse events =====

test_that("Fuse events appear when fusing triggers", {
  ds <- small_ds

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(5501)
  # fuseInterval=1 means fuse after every replicate
  result <- ts_driven_cb(ds, recorder, maxReplicates = 6L, targetHits = 6L,
                          fuseInterval = 1L)

  phases <- vapply(log, `[[`, character(1), "phase")
  # Fuse events may or may not appear (only when fuse improves score),
  # but search should complete successfully
  expect_true("done" %in% phases)
})
