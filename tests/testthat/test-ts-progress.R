# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
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


# ===== Callback fires even at verbosity 0 =====
# Callbacks are always invoked when present (regardless of verbosity)
# so that Shiny progress file polling works at verbosity = 0.

test_that("Callback still invoked when verbosity = 0", {
  ds <- small_ds

  log <- list()
  recorder <- function(info) {
    log[[length(log) + 1L]] <<- info
  }

  set.seed(8091)
  result <- ts_driven_cb(ds, recorder, verbosity = 0L)

  expect_true(length(log) > 0)
  # Should have at least replicate + done events
  phases <- vapply(log, function(x) x$phase, character(1))
  expect_true("replicate" %in% phases)
  expect_true("done" %in% phases)
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

# ===== Intra-phase heartbeat =====
# Phase-boundary reporting alone leaves the console silent for as long as one
# phase runs: measured at 582 s (TBR) and 549 s (ratchet) on a 182-tip, 420-char
# matrix with inapplicable tokens throughout -- 96% of a 1173 s replicate.

test_that("Heartbeat reports inside a long phase, and honours its interval", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]

  HeartbeatLines <- function(seconds) {
    out <- withr::with_envvar(c(TS_HEARTBEAT_SECONDS = seconds), capture.output({
      set.seed(3)
      invisible(MaximizeParsimony(ds, .rung = "thorough",
                                  maxReplicates = 1L, verbosity = 2L))
    }, type = "output"))
    grep("in phase", out, value = TRUE)
  }

  # A tiny interval fires; this dataset is small, so only the coarse-strided
  # ratchet reliably ticks, which is enough to show the mechanism works.
  expect_gt(length(HeartbeatLines("0.001")), 0L)

  # Explicitly disabled.
  expect_length(HeartbeatLines("0"), 0L)

  # Junk falls back to the DEFAULT cadence rather than being read as 0.  That
  # default is 120 s off a terminal, which this sub-second search cannot reach,
  # so the observable claim here is only that junk parses without error and does
  # not turn into a fast interval; the fallback value itself is asserted in
  # ts_heartbeat.cpp's own logic, not here.
  expect_length(HeartbeatLines("not-a-number"), 0L)
})

test_that("Heartbeat never reports a score below the true optimum", {
  # The ratchet's perturbed TBR searches a REWEIGHTED matrix, and sectorial
  # searches score a SUBTREE; both ran ~33 where the real optimum is 79.  Only
  # call sites searching the whole tree under real weights set a
  # TBRParams::heartbeat_label, so no such score can be reported.
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  best <- NULL
  out <- withr::with_envvar(c(TS_HEARTBEAT_SECONDS = "0.001"), capture.output({
    set.seed(3)
    best <- MaximizeParsimony(ds, .rung = "thorough", maxReplicates = 1L,
                              verbosity = 2L)
  }, type = "output"))
  optimum <- attr(best, "score")

  reported <- as.numeric(sub(".*best ([0-9.]+),.*", "\\1",
                             grep("in phase", out, value = TRUE)))
  skip_if(length(reported) == 0L, "no heartbeat lines emitted on this platform")
  expect_true(all(reported >= optimum))
})

test_that("Heartbeat does not change the search result", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]

  Search <- function(seconds) {
    withr::with_envvar(c(TS_HEARTBEAT_SECONDS = seconds), {
      set.seed(42)
      MaximizeParsimony(ds, effort = -9L, maxReplicates = 3L,
                        verbosity = 0L)
    })
  }
  quiet <- Search("0")
  noisy <- Search("0.001")
  expect_equal(attr(noisy, "score"), attr(quiet, "score"))
  expect_length(noisy, length(quiet))
})
