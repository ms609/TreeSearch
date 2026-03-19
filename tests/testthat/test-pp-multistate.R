# Tests for multi-state StepInformation (T-102)
context("Multi-state StepInformation")

test_that("2-state backward compatibility with LogCarter1", {
  for (a in 2:6) for (b in 2:a) {
    char <- rep(c("0", "1"), c(a, b))
    info <- StepInformation(char)
    
    # Old formula
    logProfile <- vapply(seq_len(b), LogCarter1, double(1), a, b)
    old_info <- TreeTools::Log2Unrooted(a + b) -
      (TreeSearch:::.LogCumSumExp(logProfile) / log(2))
    old_info[old_info < sqrt(.Machine[["double.eps"]])] <- 0
    names(old_info) <- seq_len(b)
    
    expect_equal(unname(info), unname(old_info),
                 tolerance = 1e-12,
                 label = paste0("(", a, ",", b, ")"))
    expect_equal(names(info), as.character(seq_len(b)),
                 label = paste0("names (", a, ",", b, ")"))
  }
})

test_that("3-state StepInformation produces correct structure", {
  # (3, 2, 2): 3 states, 7 tips, min steps = 2
  char <- rep(c("0", "1", "2"), c(3, 2, 2))
  info <- StepInformation(char)
  
  expect_true(length(info) >= 1)
  expect_equal(as.integer(names(info)[1]), 2L)  # min steps = 2
  expect_true(all(info >= 0))
  expect_true(all(diff(info) <= sqrt(.Machine[["double.eps"]])))  # monotonically decreasing
  expect_true(info[1] > 0)  # first entry has positive information
  expect_equal(unname(info[length(info)]), 0)  # last entry is zero
})

test_that("4-state StepInformation produces correct structure", {
  # (3, 2, 2, 2): 4 states, 9 tips, min steps = 3
  char <- rep(c("a", "b", "c", "d"), c(3, 2, 2, 2))
  info <- StepInformation(char)
  
  expect_true(length(info) >= 1)
  expect_equal(as.integer(names(info)[1]), 3L)
  expect_true(all(info >= 0))
  expect_true(all(diff(info) <= sqrt(.Machine[["double.eps"]])))
  expect_true(info[1] > 0)
})

test_that("infeasible multi-state falls back to binary with warning", {
  # k=5 with 11 tips: MaddisonSlatkin is infeasible; should fall back to
  # top-2 binary decomposition with a warning

  char <- rep(c("0", "1", "2", "3", "4"), c(3, 2, 2, 2, 2))
  expect_warning(info <- StepInformation(char), "reducing to 2")
  
  expect_true(length(info) >= 1)
  # Fallback uses top 2 states (3 + 2 = 5 tips), min steps = 1
  expect_equal(as.integer(names(info)[1]), 1L)
  expect_true(all(info >= 0))
  expect_true(all(diff(info) <= sqrt(.Machine[["double.eps"]])))
  expect_true(info[1] > 0)
  
  # k=3 with many tips: also falls back
  char3 <- rep(c("a", "b", "c"), c(10, 8, 6))
  expect_warning(info3 <- StepInformation(char3), "reducing to 2")
  expect_true(length(info3) >= 1)
  expect_true(all(info3 >= 0))
  
  # k=4 with many tips: also falls back
  char4 <- rep(c("x", "y", "z", "w"), c(8, 5, 4, 3))
  expect_warning(info4 <- StepInformation(char4), "reducing to 2")
  expect_true(length(info4) >= 1)
  expect_true(all(info4 >= 0))
})

test_that(">5 state warns and truncates to 5", {
  char <- rep(c("a", "b", "c", "d", "e", "f"), c(4, 3, 3, 2, 2, 2))
  # Gets two warnings: >5 truncation, then infeasibility fallback
  expect_warning(info <- StepInformation(char),
                 "5 most frequent|reducing to 2")
  
  expect_true(length(info) >= 1)
  expect_true(all(info >= 0))
  expect_true(all(diff(info) <= sqrt(.Machine[["double.eps"]])))
})

test_that("3-state with singletons includes singleton offset", {
  # 3 non-singleton states + 1 singleton state
  # split before singleton removal: c(4, 3, 2, 1), minSteps = 3
  # split after singleton removal: c(4, 3, 2), nSingletons = 1
  # reduced minSteps = 2, total minSteps = 2 + 1 = 3
  char <- rep(c("a", "b", "c", "d"), c(4, 3, 2, 1))
  info <- StepInformation(char)
  
  expect_equal(as.integer(names(info)[1]), 3L)  # 2 (reduced min) + 1 (singleton)
})

test_that("3-state information sums correctly (probabilities)", {
  # Verify that the cumulative probabilities are consistent with
  # MaddisonSlatkin summing to 1
  char <- rep(c("0", "1", "2"), c(3, 2, 2))
  info <- StepInformation(char)
  
  # At the last entry (info = 0), all trees are consistent:
  # cumsum(P) = 1, so -log2(1) = 0
  expect_equal(unname(info[length(info)]), 0)
  
  # The information should be finite and positive for early entries
  expect_true(all(is.finite(info)))
  expect_true(info[1] > 1)  # Should be substantial for a 3-state char
})

test_that("3-state matches manual MaddisonSlatkin computation", {
  # (3, 2, 2): bitmask states = c(3, 2, 0, 2, 0, 0, 0)
  states <- c(3L, 2L, 0L, 2L, 0L, 0L, 0L)
  n <- sum(states)
  logP_ms <- MaddisonSlatkin(2:(n - 1L), states)
  
  # Trim trailing -Inf
  finite <- is.finite(logP_ms)
  logP_ms <- logP_ms[seq_len(max(which(finite)))]
  
  # Cumulative info
  manual_info <- -TreeSearch:::.LogCumSumExp(logP_ms) / log(2)
  manual_info[manual_info < sqrt(.Machine[["double.eps"]])] <- 0
  
  char <- rep(c("0", "1", "2"), c(3, 2, 2))
  auto_info <- StepInformation(char)
  
  expect_equal(unname(auto_info), unname(manual_info), tolerance = 1e-12)
})

test_that("multi-state info is always >= 0", {
  test_chars <- list(
    rep(c("0", "1", "2"), c(5, 3, 2)),
    rep(c("0", "1", "2"), c(10, 5, 3)),
    rep(c("0", "1", "2", "3"), c(5, 3, 2, 2)),
    rep(c("0", "1", "2", "3", "4"), c(4, 3, 3, 2, 2))
  )
  
  for (i in seq_along(test_chars)) {
    info <- StepInformation(test_chars[[i]])
    expect_true(all(info >= 0), label = paste("test char", i))
    expect_true(all(is.finite(info)), label = paste("finite char", i))
  }
})

test_that("3-state character yields more info than binary truncation", {
  # A 3-state character should contain more information when all 3 states
  # are used than when truncated to 2. This verifies the multi-state
  # path provides additional discriminating power.
  char <- rep(c("0", "1", "2"), c(4, 3, 2))
  info_3state <- StepInformation(char)
  
  # Manual binary truncation (keep 2 largest groups)
  info_binary <- StepInformation(rep(c("0", "1"), c(4, 3)))
  
  # The multi-state character should have more total information
  # (info at minimum steps)
  expect_gt(info_3state[1], info_binary[1])
})
