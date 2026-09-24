# Implied-weights ratchet depth: constant default, coupled to the user's
# `targetHits` signal, never below the default depth, never above what was
# measured, and always yielding to an explicit `ratchetCycles`.

test_that("Implied weights deepens the ratchet under thorough and large", {
  expect_equal(
    TreeSearch:::.IwRatchetDepth("thorough", 10, 20L, 20L),
    TreeSearch:::.iwRatchetCycles
  )
  expect_equal(
    TreeSearch:::.IwRatchetDepth("large", 10, 20L, 20L),
    TreeSearch:::.iwRatchetCycles
  )
})

test_that("Equal weights and profile parsimony are left alone", {
  # Equal weights measurably does not want the deeper ratchet.
  expect_null(TreeSearch:::.IwRatchetDepth("thorough", Inf, 20L, 20L))
  # `concavity` may still be the "profile" sentinel at the call site.
  expect_null(TreeSearch:::.IwRatchetDepth("thorough", "profile", 20L, 20L))
})

test_that("Presets other than thorough and large are left alone", {
  for (preset in c("default", "sprint", "none")) {
    expect_null(TreeSearch:::.IwRatchetDepth(preset, 10, 20L, 20L))
  }
})

test_that("Raising targetHits deepens the ratchet in proportion", {
  base <- TreeSearch:::.iwRatchetCycles
  expect_equal(TreeSearch:::.IwRatchetDepth("thorough", 10, 40L, 20L), 2L * base)
  # Escalation is capped at the deepest ratchet with supporting measurements.
  expect_equal(
    TreeSearch:::.IwRatchetDepth("thorough", 10, 2000L, 20L),
    TreeSearch:::.iwRatchetMaxCycles
  )
})

test_that("Lowering targetHits does not shallow the ratchet", {
  # Shallower ratchets were slower to the optimum on every matrix tested, so the
  # documented `targetHits = 4` idiom must not drag the depth down with it.
  expect_equal(
    TreeSearch:::.IwRatchetDepth("thorough", 10, 4L, 20L),
    TreeSearch:::.iwRatchetCycles
  )
})

test_that("An explicit ratchetCycles always wins", {
  expect_null(
    TreeSearch:::.IwRatchetDepth("thorough", 10, 40L, 20L,
                                 userSet = "ratchetCycles")
  )
  # An unrelated explicit field does not suppress the coupling.
  expect_equal(
    TreeSearch:::.IwRatchetDepth("thorough", 10, 20L, 20L,
                                 userSet = c("tbrMaxHits", "driftCycles")),
    TreeSearch:::.iwRatchetCycles
  )
})

test_that("Degenerate hit counts fall back to the default depth", {
  expect_equal(
    TreeSearch:::.IwRatchetDepth("thorough", 10, 20L, 0L),
    TreeSearch:::.iwRatchetCycles
  )
  expect_equal(
    TreeSearch:::.IwRatchetDepth("thorough", 10, NA_integer_, 20L),
    TreeSearch:::.iwRatchetCycles
  )
})

test_that("A thorough implied-weights search honours ratchetCycles overrides", {
  # End-to-end: the coupling must not break the documented preset + override
  # idiom, and must leave equal-weights searches scoring as before.
  dataset <- TreeTools::MatrixToPhyDat(rbind(
    a = c(0, 0, 0, 0), b = c(0, 0, 0, 1), c = c(1, 1, 0, 0),
    d = c(1, 1, 1, 0), e = c(1, 0, 1, 1), f = c(0, 1, 1, 1)
  ))
  iw <- MaximizeParsimony(dataset, concavity = 10, .rung = "thorough",
                          maxReplicates = 2L, verbosity = 0L)
  expect_s3_class(iw, "multiPhylo")
  expect_true(all(is.finite(attr(iw, "score"))))

  pinned <- MaximizeParsimony(dataset, concavity = 10, .rung = "thorough",
                              ratchetCycles = 2L, maxReplicates = 2L,
                              verbosity = 0L)
  expect_s3_class(pinned, "multiPhylo")
  # Same optimum either way on a dataset this small: the depth changes effort,
  # never the objective.
  expect_equal(min(attr(pinned, "score")), min(attr(iw, "score")))
})
