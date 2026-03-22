## Tests for SearchControl() and the control parameter in MaximizeParsimony()

library("TreeTools")

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]

test_that("SearchControl() returns correct class and structure", {
  ctrl <- SearchControl()
  expect_s3_class(ctrl, "SearchControl")
  expect_true(is.list(ctrl))
  expect_equal(ctrl$ratchetCycles, 12L)
  expect_equal(ctrl$driftCycles, 2L)
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
  r_strict <- MaximizeParsimony(
    ds, maxReplicates = 3L, targetHits = 1L, verbosity = 0L,
    control = SearchControl(poolSuboptimal = 0)
  )
  r_sub <- MaximizeParsimony(
    ds, maxReplicates = 3L, targetHits = 1L, verbosity = 0L,
    control = SearchControl(poolSuboptimal = 5)
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
