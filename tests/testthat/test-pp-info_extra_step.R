context("pp_info_extra_step.R")

test_that("Information content of steps calculated correctly", {
  expect_equal(c(1, 2),
               as.double(ICS(2, 2, 10000, warn = FALSE) * NUnrooted(4)))
  expect_equal(c(3, 12), tolerance = 1e-5,
               as.double(ICS(2, 3, 10000, warn = FALSE) * NUnrooted(5)))
  
  
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(1)
  expect_equal(cumsum(as.double(ICS(3,3, 6000, warn = FALSE) * NUnrooted(6))),
               c(NUnrootedMult(c(3,3)),
                 NUnrootedMult(c(3,3)) + WithOneExtraStep(c(3,3)),
                 NUnrooted(6)),
               tolerance = 1)
    
  expect_equal(cumsum(as.double(ICS(3, 12, 60000, warn = FALSE))),
               c(NUnrootedMult(c(3,12)) / NUnrooted(3 + 12),
                 (NUnrootedMult(c(3,12)) + WithOneExtraStep(c(3,12)))
                 / NUnrooted(3 + 12),
                 1),
               tolerance = 5e-03)
})

test_that("WithOneExtraStep() input format", {
  expect_equal(WithOneExtraStep(7, 5), WithOneExtraStep(c(5, 7)))
})

test_that("WithOneExtraStep()", {
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
  
test_that("ExtraStep()", {
  data("profiles", package = "TreeSearch")
  
  expect_equal(c(5, 2, 0, 0, 0),
               .ExtraSteps(2, 1, 3, 1, 0, 1))
  
  expect_equal(c(7, 0, 0, 0, 0),
               .ExtraSteps(2, 1, 3, 0, 1, 2))
  
  expect_equal(c(15, 20, 0, 0, 0),
               .ExtraSteps(2, 2, 2, 1, 0, 1))
  
  expect_equal(c(15, 90, 0, 0, 0),
               .ExtraSteps(2, 3, 1, 1, 0, 1))
  
  expect_equal(c(1, 6, 0, 0, 0),
               .ExtraSteps(4, 1, 1, 1, 0, 1))
  
  expect_equal(c(0, 5, 30, 0, 0),
               .ExtraSteps(4, 2, 0, 0, 0, 0))
  
  Test <- function (a, b) {
    Profile <- function (a, b) {
      n <- sum(a, b)
      2 ^ (profiles[[n]][[2]][[n - max(a, b) - 1]] + Log2Unrooted(n))
    }
    expect_equivalent(Profile(a, b),
                      cumsum(ExtraSteps(a, b)))
  }
})