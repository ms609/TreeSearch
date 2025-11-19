test_that("MaddisonSlatkin is numerically correct", {
  expect_equal(MaddisonSlatkin(0, c(1, 1)), log(0))
  
  expect_equal(MaddisonSlatkin(1, c(1, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(2, 0)), log(1))
  expect_equal(MaddisonSlatkin(1, c(1, 0, 0, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(0, 0, 0, 2)), log(1))
})
