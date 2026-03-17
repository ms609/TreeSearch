library("TreeTools", quietly = TRUE)

make_ts_data <- function(dataset) {
  contrast <- attr(dataset, "contrast")
  tip_data <- t(vapply(dataset, as.integer, integer(length(dataset[[1]]))))
  weight <- attr(dataset, "weight")
  levels <- attr(dataset, "levels")
  list(contrast = contrast, tip_data = tip_data,
       weight = weight, levels = levels)
}

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]
td <- make_ts_data(ds)

test_that("driven search returns timings vector", {
  set.seed(4781)
  result <- TreeSearch:::ts_driven_search(
    td$contrast, td$tip_data, td$weight, td$levels,
    maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )

  expect_true("timings" %in% names(result))
  timings <- result$timings

  expected_names <- c("wagner_ms", "tbr_ms", "xss_ms", "rss_ms", "css_ms",
                      "ratchet_ms", "drift_ms", "final_tbr_ms", "fuse_ms")
  expect_equal(sort(names(timings)), sort(expected_names))
  expect_true(is.numeric(timings))

  # All timings should be non-negative
  expect_true(all(timings >= 0))

  # Core phases should have positive time
  expect_gt(timings[["wagner_ms"]], 0)
  expect_gt(timings[["tbr_ms"]], 0)
  expect_gt(timings[["ratchet_ms"]], 0)
  expect_gt(timings[["final_tbr_ms"]], 0)
})

test_that("timings sum is plausible", {
  set.seed(2198)
  result <- TreeSearch:::ts_driven_search(
    td$contrast, td$tip_data, td$weight, td$levels,
    maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )

  total <- sum(result$timings)
  expect_gt(total, 0)
  expect_lt(total, 60000)
})

test_that("MaximizeParsimony returns timings attribute", {
  set.seed(5532)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                              verbosity = 0L)

  timings <- attr(result, "timings")
  expect_false(is.null(timings))
  expect_true(is.numeric(timings))
  expect_true(all(timings >= 0))
  expect_true("wagner_ms" %in% names(timings))
  expect_true("ratchet_ms" %in% names(timings))
})

test_that("zero replicates returns zero timings", {
  result <- TreeSearch:::ts_driven_search(
    td$contrast, td$tip_data, td$weight, td$levels,
    maxReplicates = 0L, targetHits = 1L
  )

  expect_true(all(result$timings == 0))
})

test_that("timings accumulate across replicates", {
  set.seed(6401)
  result1 <- TreeSearch:::ts_driven_search(
    td$contrast, td$tip_data, td$weight, td$levels,
    maxReplicates = 1L, targetHits = 1L, verbosity = 0L
  )
  set.seed(6401)
  result3 <- TreeSearch:::ts_driven_search(
    td$contrast, td$tip_data, td$weight, td$levels,
    maxReplicates = 3L, targetHits = 3L, verbosity = 0L
  )

  # 3 replicates should have more total time than 1
  expect_gt(sum(result3$timings), sum(result1$timings) * 0.8)
})
