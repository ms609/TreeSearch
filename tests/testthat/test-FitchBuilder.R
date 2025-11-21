test_that("Fitch builder counts correctly", {
  expect_equal(AssignLeavesToRegions(c(1, 1), c(2, 2)),
               NRooted(2) * NRooted(2))
  expect_equal(AssignLeavesToRegions(c(1, 1, 1), c(3, 3, 3)),
               NRooted(3) * NRooted(3) * NRooted(3))
  
  expect_equal(BuildCounter(1, c(2, 2)), Carter1(1, 2, 2))
  
  expect_equal(BuildCounter(1, c(2, 3)), Carter1(1, 2, 3))
  expect_equal(BuildCounter(2, c(2, 3)), Carter1(2, 2, 3))
  
  expect_equal(BuildCounter(1, c(3, 3)), Carter1(1, 3, 3))
  expect_equal(BuildCounter(2, c(3, 3)), Carter1(2, 3, 3))
  expect_equal(BuildCounter(3, c(3, 3)), Carter1(3, 3, 3))
  
  
  
})
