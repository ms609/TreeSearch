# Fractional per-character weights via attr(dat, "weight").
#
# The C++ scoring engine still stores int weights, so .ScaleWeight() converts
# fractional vectors to integer at the R-level chokepoints with a configurable
# precision (option TreeSearch.fractional.scale, default 1260 = 2*2*3*3*5*7).

weight_multiplier <- 1260L

test_that(".ScaleWeight is a no-op for integer weights", {
  w <- c(1L, 2L, 3L, 1L)
  expect_identical(TreeSearch:::.ScaleWeight(w), w)
})

test_that(".ScaleWeight passes integer-valued doubles through unscaled", {
  w <- c(1, 2, 3)
  expect_identical(TreeSearch:::.ScaleWeight(w), c(1L, 2L, 3L))
})

test_that(".ScaleWeight scales true fractional weights by default", {
  w <- c(0.5, 1.25, 2.0)
  expect_identical(TreeSearch:::.ScaleWeight(w),
                   as.integer(round(w * weight_multiplier)))
})

test_that(".ScaleWeight honours TreeSearch.fractional.scale option", {
  old <- options(TreeSearch.fractional.scale = 100L)
  on.exit(options(old), add = TRUE)
  expect_identical(TreeSearch:::.ScaleWeight(c(0.5, 0.75)),
                   c(50L, 75L))
})

test_that(".ScaleWeight floors tiny positive weights at 1 to avoid drop", {
  # A weight whose scaled value rounds to 0 would silently drop the character.
  # Floor at 1 to preserve it.
  tiny <- 0.4 / weight_multiplier  # scales to 0 without the floor
  expect_identical(TreeSearch:::.ScaleWeight(c(tiny, 0.5)),
                   c(1L, as.integer(round(0.5 * weight_multiplier))))
})

test_that("TreeLength respects fractional weights end-to-end", {
  dat <- TreeTools::MatrixToPhyDat(matrix(
    c(0, 0, 1, 1,
      0, 1, 0, 1,
      0, 1, 1, 0,
      1, 0, 0, 1),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"),
                    c("c1", "c2", "c3", "c4"))
  ))
  tr <- ape::read.tree(text = "((A, B), (C, D));")

  base <- TreeLength(tr, dat)  # integer weights, all 1

  attr(dat, "weight") <- c(2, 2, 2, 2)  # numeric but integer-valued
  expect_equal(TreeLength(tr, dat), base * 2)

  attr(dat, "weight") <- rep(0.5, 4L)
  half <- TreeLength(tr, dat)
  expect_equal(half, base * 0.5 * weight_multiplier)

  attr(dat, "weight") <- rep(1.5, 4L)
  one_half <- TreeLength(tr, dat)
  expect_equal(one_half, base * 1.5 * weight_multiplier)
})

test_that(".ScaleWeight errors when sum(scaled) > .Machine$integer.max", {
  # Each weight of (INT_MAX / 4 + 1) * scale would push total >> INT_MAX.
  # Use a non-integer value so the fractional branch runs.
  big_w <- (.Machine$integer.max %/% 4L + 1L) / weight_multiplier
  old <- options(TreeSearch.fractional.scale = weight_multiplier)
  on.exit(options(old), add = TRUE)
  expect_error(
    TreeSearch:::.ScaleWeight(rep(big_w, 5L)),
    regexp = "integer.max",
    fixed = FALSE
  )
})

test_that("Resample() errors cleanly when sum(weights) > INT_MAX", {
  # Bypass .ScaleWeight() by setting integer weights directly on the phyDat.
  # This exercises the C++ guard in ts_resample.cpp.
  big_w <- .Machine$integer.max %/% 4L + 1L
  dat <- TreeTools::MatrixToPhyDat(matrix(
    c(0, 0, 1, 1, 0,
      0, 1, 0, 1, 0,
      0, 1, 1, 0, 1,
      1, 0, 0, 1, 0,
      1, 0, 1, 0, 0),
    nrow = 5, byrow = TRUE,
    dimnames = list(letters[1:5], paste0("c", 1:5))
  ))
  # Force integer weights summing to > INT_MAX on the nr unique patterns.
  attr(dat, "weight") <- rep(big_w, attr(dat, "nr"))
  expect_error(
    Resample(dat, nReplicates = 1L),
    regexp = "INT_MAX|integer.max",
    ignore.case = TRUE
  )
})

test_that("MaximizeParsimony accepts fractional weights without crashing", {
  set.seed(42)
  dat <- TreeTools::MatrixToPhyDat(matrix(
    sample(0:1, 60, replace = TRUE),
    nrow = 6, byrow = TRUE,
    dimnames = list(letters[1:6], paste0("c", 1:10))
  ))
  attr(dat, "weight") <- runif(attr(dat, "nr"), min = 0.1, max = 1)
  expect_error(
    MaximizeParsimony(dat, maxSeconds = 5L, maxReplicates = 4L,
                      verbosity = 0L),
    NA
  )
})
