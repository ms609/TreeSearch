# The no-improvement stopping rule computes (targetHits / hits) * nTip * perturbStopFactor.
# Large-but-legal settings overflow `int`; casting the out-of-range double was undefined
# behaviour that landed negative, so the rule fired on the FIRST non-improving replicate.
# The symptom was silent: the search returned a worse tree, with no warning, from settings
# that ask for MORE search rather than less.

test_that("A huge targetHits x perturbStopFactor product does not stop the search", {
  dataset <- TreeTools::MatrixToPhyDat(rbind(
    a = c(0, 0, 0, 0, 1, 0), b = c(0, 0, 0, 1, 1, 1), c = c(1, 1, 0, 0, 0, 1),
    d = c(1, 1, 1, 0, 1, 0), e = c(1, 0, 1, 1, 0, 0), f = c(0, 1, 1, 1, 0, 1),
    g = c(1, 0, 0, 1, 1, 0), h = c(0, 1, 0, 0, 1, 1)
  ))
  reps <- 12L
  # 99999 * 8 tips * 10000 overflows int; a negative limit stops after ~2 replicates.
  r <- MaximizeParsimony(dataset, targetHits = 99999L, perturbStopFactor = 10000L,
                         maxReplicates = reps, verbosity = 0L)
  expect_s3_class(r, "multiPhylo")
  # Neither stopping rule can fire here (99999 hits is unreachable in 12 replicates), so
  # the replicate cap is what must end the search.
  expect_equal(attr(r, "replicates"), reps)
})

test_that("The no-improvement rule still fires when it genuinely should", {
  # Guards the saturation against being a blanket disable: a small factor must still stop
  # the search well short of its replicate cap.
  dataset <- TreeTools::MatrixToPhyDat(rbind(
    a = c(0, 0, 0, 0), b = c(0, 0, 0, 1), c = c(1, 1, 0, 0),
    d = c(1, 1, 1, 0), e = c(1, 0, 1, 1), f = c(0, 1, 1, 1)
  ))
  r <- MaximizeParsimony(dataset, targetHits = 2L, perturbStopFactor = 1L,
                         maxReplicates = 500L, verbosity = 0L)
  expect_s3_class(r, "multiPhylo")
  expect_lt(attr(r, "replicates"), 500L)
})
