test_that("Downpass handles overlap vs union correctly", {
  dp <- .DownpassOutcome(4)
  
  # {A}, {A,B} -> {A}
  expect_equal(dp[1, 3], as.raw(1)) 
  
  # {A}, {B} -> {A,B}
  expect_equal(dp[1, 2], as.raw(3)) 
  
  # {A,B}, {A,B,C} -> {A,B}
  expect_equal(dp[3, 7], as.raw(3)) 
  
  # symmetry
  expect_equal(dp[2, 1], dp[1, 2])
})

test_that("Uppass covers all cases", {
  up <- .UppassOutcome(4)
  
  # Case I: ancestor subset of node
  expect_equal(up[1, 3, 1], as.raw(1)) # {A}, {A,B}, ancestor {A}
  
  # Case III: disjoint descendants
  expect_equal(up[1, 2, 4], as.raw(7)) # {A}, {B}, ancestor {C}
  
  # Default case
  expect_equal(up[3, 6, 4], as.raw(6)) # {A,B}, {B,C}, ancestor {C} => {B,C}
  
  # Ancestor irrelevant
  expect_equal(up[1, 1, 8], as.raw(1)) # {A}, {A}, ancestor {D}
})
