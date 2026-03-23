context("pp_info_extra_step.R")
library("TreeSearch", quietly = TRUE)

test_that("Bad input safely handled", {
  expect_equal(0, WithOneExtraStep(1))
  expect_error(WithOneExtraStep(2, 2, 2))
  
  expect_equal(0, Carter1(5, 6, 4))
  expect_equal(-Inf, LogCarter1(5, 6, 4))
  expect_equal(-Inf, Log2Carter1(5, 6, 4))
})

test_that("StepInformation() works", {
  expect_equal(c(`0` = 0), StepInformation(rep(3L, 10), ambiguousTokens = 3))
  expect_equal(c(`0` = 0), StepInformation(c(4L, rep(3L, 10)), 3))
  expect_true(all(is.finite(StepInformation(rep.int(1:3, times = c(139, 45, 41)),
                                            ambiguousTokens = 3))))
  expect_true(all(is.finite(StepInformation(
    char = rep.int(1:2, times = c(600, 600))))))
})

test_that("Carter1() matches profile counts", {
  data("profiles", package = "TreeSearch")
  Test <- function (a, b) {
    n <- sum(a, b)
    counted <- 2 ^ profiles[[n]][[2]][[n - max(a, b) - 1]] * NUnrooted(n)
    m <- as.integer(names(counted))
    for (mi in m) {
      expect_equal(log2(Carter1(mi, a, b)), Log2Carter1(mi, b, a))
      expect_equal(log(Carter1(mi, a, b)), LogCarter1(mi, b, a))
    }
    expect_equivalent(counted,
                      cumsum(vapply(m, Carter1, a = a, b = b, double(1))))
  }
  
  Test(2, 4)
  Test(2, 5)
  Test(2, 6)
  Test(2, 7)
  Test(2, 8)
  
  Test(3, 4)
  Test(3, 5)
  Test(3, 6)
  Test(3, 7)
  
  Test(4, 4)
  Test(4, 5)
  Test(4, 6)
  
  Test(5, 4)
  Test(5, 5)
  
})

test_that("WithOneExtraStep() input format", {
  expect_equal(WithOneExtraStep(7, 5), WithOneExtraStep(c(5, 7)))
})

test_that("WithOneExtraStep()", {
  library("TreeTools", quietly = TRUE)
  data("profiles", package = "TreeSearch")
  Test <- function (a, b) {
    n <- sum(a, b)
    expect_equivalent(2 ^ profiles[[n]][[2]][[n - max(a, b) - 1]][2] * NUnrooted(n),
                      NUnrootedMult(c(a, b)) + WithOneExtraStep(c(a, b)))
  }
  
  Test(4, 2)
  Test(3, 3)
  Test(8, 2)
  Test(4, 3)
  Test(7, 3)
  Test(6, 4)
  Test(5, 5)
  
  expect_equal(NUnrooted(6) / NUnrooted(5) * WithOneExtraStep(2:3),
               WithOneExtraStep(1:3))
  expect_equal(NUnrooted(10) / NUnrooted(5) * WithOneExtraStep(2:3),
               WithOneExtraStep(2:3, rep(1, 5)))
})

test_that(".LogCumSumExp()", {
  Test <- function (x) {
    naive <- log(cumsum(exp(x)))
    if (all(is.finite(naive))) {
      expect_equal(naive, TreeSearch:::.LogCumSumExp(x))
    } else {
      expect_true(all(is.finite(TreeSearch:::.LogCumSumExp(x))))
    }
  }
  Test(log(c(1:5, 5:1)))
  Test(c(10, 700, 100))
  Test(c(10, 7000, 100))
})

test_that(".LogCumSumExp() handles -Inf without NaN", {
  # Both x[k] and Lk[k-1] = -Inf: IEEE 754 gives -Inf - (-Inf) = NaN.
  # Guard must keep result = -Inf (log(0 + 0) = -Inf), not NaN.
  out <- TreeSearch:::.LogCumSumExp(c(-Inf, -Inf, -3.0))
  expect_false(any(is.nan(out)), label = "no NaN when consecutive -Inf")
  expect_equal(out[1], -Inf)
  expect_equal(out[2], -Inf)
  expect_equal(out[3], -3.0)   # log(exp(-3)) = -3

  # Single -Inf at start, then finite: should recover correctly
  out2 <- TreeSearch:::.LogCumSumExp(c(-Inf, -2.0, -5.0))
  expect_false(any(is.nan(out2)))
  expect_equal(out2[1], -Inf)
  expect_equal(out2[2], -2.0)

  # All -Inf: result should be all -Inf, not NaN
  out3 <- TreeSearch:::.LogCumSumExp(c(-Inf, -Inf, -Inf))
  expect_true(all(out3 == -Inf))
  expect_false(any(is.nan(out3)))
})