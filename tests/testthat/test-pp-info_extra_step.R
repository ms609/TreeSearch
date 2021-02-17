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

test_that("Carter1() matches profile counts", {
  data("profiles", package = "TreeSearch")
  Test <- function (a, b) {
    n <- sum(a, b)
    counted <- 2 ^ profiles[[n]][[2]][[n - max(a, b) - 1]] * NUnrooted(n)
    m <- as.integer(names(counted))
    for (mi in m) {
      expect_equal(log2(Carter1(mi, a, b)), Log2Carter1(mi, b, a))
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
  
