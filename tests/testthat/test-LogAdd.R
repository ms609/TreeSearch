test_that("LogAdd works", {
  expect_log_add <- function(x, y) {
    expect_equal(LogAdd(log(x), log(y)), log(x + y))
  }
  expect_log_add(0, 1)
  expect_log_add(0, 0)
  expect_log_add(1, 0)
  expect_log_add(1, 1)
  expect_log_add(c(0, 0, 1, 1), c(1, 0, 0, 1))
  expect_log_add(1:4, 4:1)
  expect_log_add(1e-9, 1e-9)
  expect_equal(LogAdd(log(1:3), log(0:2), log(5:3)), log(1:3 + 0:2 + 5:3))
})
