# Fractional per-character weights via attr(dat, "weight").
#
# The C++ scoring engine still stores int weights, so .ScaleWeight() converts
# fractional vectors to integer at the R-level chokepoints with a configurable
# precision (option TreeSearch.fractional.scale, default 1000).

test_that(".ScaleWeight is a no-op for integer weights", {
  w <- c(1L, 2L, 3L, 1L)
  expect_identical(TreeSearch:::.ScaleWeight(w), w)
})

test_that(".ScaleWeight passes integer-valued doubles through unscaled", {
  w <- c(1, 2, 3)
  expect_identical(TreeSearch:::.ScaleWeight(w), c(1L, 2L, 3L))
})

test_that(".ScaleWeight scales true fractional weights by 1000", {
  w <- c(0.5, 1.25, 2.0)
  expect_identical(TreeSearch:::.ScaleWeight(w),
                   c(500L, 1250L, 2000L))
})

test_that(".ScaleWeight honours TreeSearch.fractional.scale option", {
  withr::with_options(list(TreeSearch.fractional.scale = 100L), {
    expect_identical(TreeSearch:::.ScaleWeight(c(0.5, 0.75)),
                     c(50L, 75L))
  })
})

test_that(".ScaleWeight floors tiny positive weights at 1 to avoid drop", {
  # 0.0005 * 1000 rounds to 0; without floor the character would silently
  # vanish from scoring. Floor to 1.
  expect_identical(TreeSearch:::.ScaleWeight(c(0.0005, 0.5)),
                   c(1L, 500L))
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

  attr(dat, "weight") <- rep(0.5, 4L)  # fractional, scale=1000 -> 500 each
  scaled_500 <- TreeLength(tr, dat)
  expect_equal(scaled_500, base * 500)  # half-weight x scale=1000 = 500 x

  attr(dat, "weight") <- rep(1.5, 4L)
  scaled_1500 <- TreeLength(tr, dat)
  expect_equal(scaled_1500, base * 1500)
})

test_that("MaximizeParsimony accepts fractional weights without crashing", {
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
