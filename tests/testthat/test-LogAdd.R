test_that("LogAdd works", {
  expect_log_add <- function(...) {
    expect_equal(LogSumExp(log(c(...))), log(sum(...)))
  }
  expect_log_add(0, 1)
  expect_log_add(0, 0)
  expect_log_add(1, 0)
  expect_log_add(1, 1)
  expect_log_add(c(0, 0, 1, 1), c(1, 0, 0, 1))
  expect_log_add(1:4, 4:1)
  expect_log_add(1e-9, 1e-9)
})
