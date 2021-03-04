test_that("Maddison Slatkin implemented successfully", {
  
  undefined <- 0
  
  expect_equal(2/3, .B0(3, 1))
  expect_equal(1/3, .B01(3, 1))
  expect_equal(0, .B1(3, 1))
  
  expect_equal(0, .P0(0, 3, 1))
  expect_equal(0, .P1(0, 3, 1))
  expect_equal(0, .P01(0, 3, 1))
  
  expect_equal(1, .P0(1, 3, 1))
  expect_equal(undefined, .P1(1, 3, 1))
  expect_equal(1, .P01(1, 3, 1))
  
  
  expect_equal(0, .P0(0, 3, 2))
  expect_equal(0, .P1(0, 3, 2))
  expect_equal(0, .P01(0, 3, 2))
  expect_equal(undefined, .P0(1, 3, 2))
  expect_equal(1, .P1(1, 3, 2))
  expect_equal(1, .P01(1, 3, 2))
  
  
  TestBSum <- function (n, i) {
    expect_equal(1, sum(.B0(n, i), .B01(n, i), .B1(n, i)))
  }
  TestPSum <- function (n, i) {
    expect_equal(1, sum(vapply(0:i, function (s) 
      sum(.P0(s, n, i) * .B0(n, i),
          .P01(s, n, i) * .B01(n, i),
          .P1(s, n, i) * .B1(n, i)), double(1))))
  }
  TestBSum(3, 0)
  TestBSum(3, 1)
  TestBSum(3, 2)
  TestBSum(3, 3)
  
  TestPSum(3, 0)
  TestPSum(3, 1)
  TestPSum(3, 2)
  TestPSum(3, 3)
  
  
  for (i in 0:8) TestBSum(8, i)
  for (i in 0:8) TestPSum(8, i)
  
  expect_equal(MaddisonSlatkin(2, 4) * NUnrooted(6),
               vapply(1:2, Carter1, double(1), 2, 4))
  
  expect_equal(MaddisonSlatkin(8, 8) * NUnrooted(16),
               vapply(1:8, Carter1, double(1), 8, 8))
  
  expect_equal(MaddisonSlatkin(12, 4) * NUnrooted(16),
               vapply(1:4, Carter1, double(1), 12, 4))
  
})
