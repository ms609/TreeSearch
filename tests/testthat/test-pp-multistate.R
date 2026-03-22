# Tests for multi-state StepInformation
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

test_that("infeasible multi-state uses MC preserving all states", {
  # Feasibility uses partition-aware split_count (sc).
  # Thresholds: k=3 sc>75, k=4 sc>50, k=5 sc>35.
  # Infeasible characters now use MC approximation (no state reduction).

  set.seed(6391)

  # k=3 n=38 (13,13,12): sc=140 >> 75
  char3 <- rep(c("a", "b", "c"), c(13, 13, 12))
  info3 <- StepInformation(char3, n_mc = 5000L)
  expect_true(length(info3) >= 1)
  # MC preserves 3 states: min steps = k - 1 = 2

  expect_equal(as.integer(names(info3)[1L]), 2L)
  expect_true(all(info3 >= 0))
  expect_true(all(diff(info3) <= sqrt(.Machine[["double.eps"]])))

  # k=4 n=24 (7,6,6,5): sc=224 >> 50
  char4 <- rep(c("x", "y", "z", "w"), c(7, 6, 6, 5))
  info4 <- StepInformation(char4, n_mc = 5000L)
  expect_true(length(info4) >= 1)
  expect_equal(as.integer(names(info4)[1L]), 3L)
  expect_true(all(info4 >= 0))

  # k=5 n=15 (4,3,3,3,2): sc=143 >> 35
  char5 <- rep(c("0", "1", "2", "3", "4"), c(4, 3, 3, 3, 2))
  info5 <- StepInformation(char5, n_mc = 5000L)
  expect_true(length(info5) >= 1)
  expect_equal(as.integer(names(info5)[1L]), 4L)
  expect_true(all(info5 >= 0))
})

test_that("approx='mc' matches exact within 1 bit for feasible character", {
  # k=3 n=7: exact is fast; compare MC approximation to exact
  char <- rep(c("0", "1", "2"), c(3, 2, 2))
  
  set.seed(4412)
  info_exact <- StepInformation(char, approx = "exact")
  info_mc    <- StepInformation(char, approx = "mc", n_mc = 10000L)
  
  common <- intersect(names(info_exact), names(info_mc))
  expect_true(length(common) >= 1L)
  # MC should agree with exact within 1 bit at every step count
  expect_true(all(abs(info_mc[common] - info_exact[common]) <= 1),
              label = "MC agrees with exact within 1 bit")
  
  # Both are non-negative and non-increasing
  expect_true(all(info_mc >= 0))
  expect_true(all(diff(info_mc) <= sqrt(.Machine[["double.eps"]])))
})

test_that("approx='mc' returns multi-state step range for infeasible char", {
  # k=3 n=38 (13,13,12): sc=140 >> 75 threshold, infeasible for exact.
  # MC should return IC starting at step 2 (k-1), not step 1 (binary)
  char <- rep(c("0", "1", "2"), c(13, 13, 12))
  
  set.seed(7731)
  info_mc <- StepInformation(char, approx = "mc", n_mc = 2000L)
  
  expect_true(length(info_mc) >= 1L)
  # Min steps = k - 1 = 2 (not 1 as binary fallback would give)
  expect_equal(as.integer(names(info_mc)[1L]), 2L)
  expect_true(all(info_mc >= 0))
  expect_true(all(diff(info_mc) <= sqrt(.Machine[["double.eps"]])))
  expect_true(info_mc[1L] > 0)
})

test_that("PrepareDataProfile preserves multi-state patterns", {
  # A small dataset with one 3-state char with many tips
  set.seed(3058)
  n <- 20L
  nchar <- 5L
  mat <- matrix(
    c(rep(0:2, c(8L, 7L, 5L)),           # char 1: 3-state, n=20 (feasible; sc=42)
      sample(0:1, n * (nchar - 1L), replace = TRUE)),
    nrow = n,
    dimnames = list(paste0("t", seq_len(n)), paste0("c", seq_len(nchar)))
  )
  dat <- TreeTools::MatrixToPhyDat(mat)

  # "auto" and "mc" should produce identical structure (both use MC for
  # infeasible chars, exact for feasible ones)
  info_auto <- PrepareDataProfile(dat, approx = "auto", n_mc = 5000L)
  info_mc   <- PrepareDataProfile(dat, approx = "mc", n_mc = 5000L)

  auto_steps <- nrow(attr(info_auto, "info.amounts"))
  mc_steps   <- nrow(attr(info_mc,   "info.amounts"))
  expect_equal(auto_steps, mc_steps)

  # Both produce valid, finite info.amounts
  expect_true(all(is.finite(attr(info_auto, "info.amounts"))))
  expect_true(all(is.finite(attr(info_mc,   "info.amounts"))))
})

test_that(">5 state warns and truncates to 5", {
  char <- rep(c("a", "b", "c", "d", "e", "f"), c(4, 3, 3, 2, 2, 2))
  expect_warning(info <- StepInformation(char, n_mc = 5000L),
                 "5 most frequent")

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
  set.seed(8203)
  test_chars <- list(
    rep(c("0", "1", "2"), c(5, 3, 2)),
    rep(c("0", "1", "2"), c(10, 5, 3)),
    rep(c("0", "1", "2", "3"), c(5, 3, 2, 2)),
    rep(c("0", "1", "2", "3", "4"), c(4, 3, 3, 2, 2))
  )

  for (i in seq_along(test_chars)) {
    info <- suppressWarnings(StepInformation(test_chars[[i]], n_mc = 5000L))
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

test_that("MC approximation matches exact within 2 bits at boundary", {
  # (5,5,5) n=15, k=3: feasible (sc=27 < 75), so exact is available.
  # Compare MC to exact to validate the log-quadratic interpolation.
  char <- rep(c("0", "1", "2"), c(5, 5, 5))

  info_exact <- StepInformation(char, approx = "exact")

  set.seed(5072)
  info_mc <- StepInformation(char, approx = "mc", n_mc = 50000L)

  common <- intersect(names(info_exact), names(info_mc))
  expect_true(length(common) >= 1L)

  # MC should agree with exact within 2 bits at every step count
  diffs <- abs(info_mc[common] - info_exact[common])
  expect_true(all(diffs <= 2),
              label = paste("max MC-exact diff:", round(max(diffs), 3), "bits"))

  # IC(0) should be close: exact P(s_min) is the same in both
  expect_equal(unname(info_mc[1L]), unname(info_exact[1L]),
               tolerance = 0.5)
})

test_that("log-quadratic interpolation produces monotone IC", {
  # Infeasible 3-state: (13,13,12), n=38. MC must produce monotonically
  # decreasing IC (non-increasing).
  set.seed(2849)
  char <- rep(c("0", "1", "2"), c(13, 13, 12))
  info <- StepInformation(char, n_mc = 10000L)

  expect_true(length(info) >= 1L)
  expect_true(all(info >= 0))
  expect_true(all(diff(info) <= sqrt(.Machine[["double.eps"]])))
  # First entry has positive information
  expect_true(info[1L] > 0)
})
