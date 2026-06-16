## Tests for SearchControl() and the control parameter in MaximizeParsimony()

library("TreeTools")

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]

test_that("SearchControl() returns correct class and structure", {
  ctrl <- SearchControl()
  expect_s3_class(ctrl, "SearchControl")
  expect_true(is.list(ctrl))
  expect_equal(ctrl$ratchetCycles, 12L)
  expect_equal(ctrl$driftCycles, 0L)
  expect_equal(ctrl$poolSuboptimal, 0)
  expect_false(ctrl$sprFirst)
  expect_equal(ctrl$ratchetPerturbProb, 0.25)
  expect_equal(ctrl$ratchetPerturbMaxMoves, 5L)
  expect_equal(ctrl$driftAfdLimit, 5L)
  expect_equal(ctrl$driftRfdLimit, 0.15)
})

test_that("SearchControl() accepts custom values", {
  ctrl <- SearchControl(ratchetCycles = 20L, driftCycles = 0L,
                         poolSuboptimal = 3.5)
  expect_equal(ctrl$ratchetCycles, 20L)
  expect_equal(ctrl$driftCycles, 0L)
  expect_equal(ctrl$poolSuboptimal, 3.5)
  # Other fields keep defaults
 expect_equal(ctrl$xssRounds, 3L)
})

test_that("SearchControl() rejects crash-inducing count parameters (RT-CPP-01)", {
  # Zero partitions -> integer division by zero in xss_partition() (SIGFPE);
  # zero poolMaxSize -> out-of-bounds read in TreePool::add() (segfault). Both
  # are uncatchable crashes, so they must be rejected at the R boundary.
  expect_error(SearchControl(xssPartitions = 0L),  "xssPartitions.*positive")
  expect_error(SearchControl(xssPartitions = -1L), "xssPartitions.*positive")
  expect_error(SearchControl(xssPartitions = NA_integer_), "xssPartitions.*positive")
  expect_error(SearchControl(cssPartitions = 0L),  "cssPartitions.*positive")
  expect_error(SearchControl(poolMaxSize = 0L),    "poolMaxSize.*positive")
  expect_error(SearchControl(poolMaxSize = -3L),   "poolMaxSize.*positive")
  # The boundary value 1 (single sector / single-tree pool) is valid.
  expect_equal(SearchControl(xssPartitions = 1L)$xssPartitions, 1L)
  expect_equal(SearchControl(cssPartitions = 1L)$cssPartitions, 1L)
  expect_equal(SearchControl(poolMaxSize = 1L)$poolMaxSize, 1L)
})

test_that("search survives a non-positive poolMaxSize via a raw control list", {
  # A plain list control bypasses SearchControl()'s validation, so the C++
  # TreePool must itself clamp max_size >= 1 (else entries_[0] on an empty
  # pool segfaults). If the clamp regresses this crashes the worker, which is
  # an acceptable loud signal for so severe a bug.
  skip_on_cran()
  set.seed(1)
  dat <- TreeTools::MatrixToPhyDat(matrix(
    sample(0:1, 24 * 30, replace = TRUE), nrow = 24,
    dimnames = list(paste0("t", 1:24), NULL)))
  ctrl <- SearchControl(ratchetCycles = 1L)
  ctrl$poolMaxSize <- 0L
  res <- MaximizeParsimony(dat, maxReplicates = 2L,
                           verbosity = 0L, control = ctrl)
  expect_true(is.finite(attr(res, "score")))
})

test_that("print.SearchControl works", {
  ctrl <- SearchControl()
  expect_output(print(ctrl), "SearchControl object")
  expect_output(print(ctrl), "Ratchet:")
  expect_output(print(ctrl), "ratchetCycles")
})

test_that("MaximizeParsimony accepts control parameter", {
  set.seed(8472)
  result <- MaximizeParsimony(
    ds, maxReplicates = 2L, targetHits = 1L, verbosity = 0L,
    control = SearchControl(ratchetCycles = 2L, driftCycles = 0L)
  )
  expect_s3_class(result, "multiPhylo")
  expect_true(attr(result, "score") > 0)
})

test_that("MaximizeParsimony accepts plain list as control", {
  set.seed(8472)
  result <- MaximizeParsimony(
    ds, maxReplicates = 2L, targetHits = 1L, verbosity = 0L,
    control = list(ratchetCycles = 2L, driftCycles = 0L)
  )
  expect_s3_class(result, "multiPhylo")
  expect_true(attr(result, "score") > 0)
})

test_that("Backward compat: ... overrides control defaults", {
  set.seed(8472)
  result <- MaximizeParsimony(
    ds, maxReplicates = 2L, targetHits = 1L, verbosity = 0L,
    ratchetCycles = 2L, driftCycles = 0L
  )
  expect_s3_class(result, "multiPhylo")
  expect_true(attr(result, "score") > 0)
})

test_that("Strategy preset overrides SearchControl defaults", {
  set.seed(8472)
  # sprint preset sets driftCycles=0, ratchetCycles=3
  r1 <- MaximizeParsimony(
    ds, strategy = "sprint",
    maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )
  expect_s3_class(r1, "multiPhylo")
})

test_that("Explicit control overrides strategy preset", {
  set.seed(8472)
  # sprint preset sets ratchetCycles=3; override to 1
  r1 <- MaximizeParsimony(
    ds, strategy = "sprint",
    maxReplicates = 2L, targetHits = 1L, verbosity = 0L,
    control = SearchControl(ratchetCycles = 1L)
  )
  expect_s3_class(r1, "multiPhylo")
})

test_that("poolSuboptimal via control collects suboptimal trees", {
  set.seed(8472)
  # Use unlimited outer resets for consistent search depth across platforms
  r_strict <- MaximizeParsimony(
    ds, maxReplicates = 3L, targetHits = 1L, verbosity = 0L,
    control = SearchControl(poolSuboptimal = 0, maxOuterResets = -1L)
  )
  r_sub <- MaximizeParsimony(
    ds, maxReplicates = 3L, targetHits = 1L, verbosity = 0L,
    control = SearchControl(poolSuboptimal = 5, maxOuterResets = -1L)
  )
  expect_gte(length(r_sub), length(r_strict))
})

test_that("Unknown ... args produce warning", {
  expect_warning(
    MaximizeParsimony(ds, maxReplicates = 1L, targetHits = 1L,
                       verbosity = 0L, fakeParam = 42),
    "Unknown"
  )
})

test_that("IW via control works", {
  set.seed(8472)
  result <- MaximizeParsimony(
    ds, concavity = 10,
    maxReplicates = 2L, targetHits = 1L, verbosity = 0L,
    control = SearchControl(ratchetCycles = 2L)
  )
  score <- attr(result, "score")
  expect_true(is.finite(score) && score > 0)
  # Verify against TreeLength
  tl <- TreeLength(result[[1]], ds, concavity = 10)
  expect_equal(score, tl, tolerance = 0.01)
})

test_that("stallEscalateFactor is validated and stored", {
  expect_equal(SearchControl()$stallEscalateFactor, 1)
  expect_equal(SearchControl(stallEscalateFactor = 2.5)$stallEscalateFactor, 2.5)
  # A value < 1 would shrink perturbation on stalling (wrong direction); NA and
  # non-scalar are also rejected at the R boundary.
  expect_error(SearchControl(stallEscalateFactor = 0.5),
               "stallEscalateFactor.*>= 1")
  expect_error(SearchControl(stallEscalateFactor = NA_real_),
               "stallEscalateFactor.*>= 1")
  expect_error(SearchControl(stallEscalateFactor = c(1, 2)),
               "stallEscalateFactor.*>= 1")
})

test_that("stallEscalateFactor escalates on stall and still scores correctly", {
  skip_on_cran()
  # Enough replicates that the search stalls (no improvement for >= nTip/10
  # reps), which engages the escalator's stalled branch; the early reps exercise
  # the non-stalled branch. Escalation alters the perturbation trajectory, not
  # the achievable optimum, so the result is a valid, correctly-scored tree.
  set.seed(8472)
  result <- MaximizeParsimony(
    ds, maxReplicates = 8L, targetHits = 99L, verbosity = 0L,
    control = SearchControl(stallEscalateFactor = 3, ratchetCycles = 3L)
  )
  expect_s3_class(result, "multiPhylo")
  score <- attr(result, "score")
  expect_true(is.finite(score) && score > 0)
  expect_equal(score, TreeLength(result[[1]], ds), tolerance = 0.01)
})
