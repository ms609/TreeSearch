# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

library("TreeTools", quietly = TRUE)

test_that("MaddisonSlatkin() recursion bottoms", {
  expect_equal(MaddisonSlatkin(0, c(1, 1)), log(0))
  expect_equal(MaddisonSlatkin(1, c(1, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(2, 0)), log(1))
  expect_equal(MaddisonSlatkin(1, c(1, 0, 0, 1)), log(1))
  expect_equal(MaddisonSlatkin(0, c(0, 0, 0, 2)), log(1))
})

test_that("MaddisonSlatkin() brute-force matches small trees", {
  # Enumerate ALL unrooted trees, count actual step distribution,
  # verify MaddisonSlatkin matches exactly.
  expect_slatkin <- function(tokens) {
    ch <- rep(seq_along(tokens), tokens)
    nTaxa <- length(ch)
    phyChar <- StringToPhyDat(paste0(ch, collapse = ""))
    trees <- as.phylo(seq_len(NUnrooted(nTaxa)) - 1L, nTaxa)
    counts <- vapply(trees, TreeLength, double(1), phyChar) |>
      tabulate()
    out <- vapply(seq_along(counts), MaddisonSlatkin, double(1),
                  tabulate(ch)) |>
      exp() * length(trees)
    expect_equal(out, counts)
  }

  # 2-state cases
  expect_slatkin(c(2, 2))     # 4 tips, 3 trees
  expect_slatkin(c(2, 3))     # 5 tips, 15 trees
  expect_slatkin(c(2, 4))     # 6 tips, 105 trees

  # 3-state case (states at bitmask positions 1, 2, 4)
  expect_slatkin(c(2, 3, 0, 2))  # 7 tips, 945 trees
})

test_that("MaddisonSlatkin() matches published examples", {
  # Maddison & Slatkin (1991) Table 1 cross-validation
  expect_equal(MaddisonSlatkin(1, c(8, 24)) + LnUnrooted(32),
               LogCarter1(1, 8, 24))
  expect_equal(MaddisonSlatkin(2, c(8, 24)) + LnUnrooted(32),
               LogCarter1(2, 8, 24))
  expect_equal(MaddisonSlatkin(3, c(7, 18)) + LnUnrooted(25),
               LogCarter1(3, 7, 18))
})

test_that("MaddisonSlatkin matches LogCarter1 for 2-state", {
  # Character: 3 leaves with state 0, 3 with state 1
  states <- c(3L, 3L, 0L)
  ms <- MaddisonSlatkin(1:3, states)
  lc <- vapply(1:3, LogCarter1, double(1), 3, 3)
  lnTotal <- TreeTools::LnUnrooted(6)

  # MaddisonSlatkin returns log(fraction); LogCarter1 returns log(count)
  expect_equal(ms + lnTotal, lc, tolerance = 1e-12)
})

test_that("MaddisonSlatkin matches LogCarter1 for asymmetric 2-state", {
  states <- c(5L, 2L, 0L)
  ms <- MaddisonSlatkin(1:2, states)
  lc <- vapply(1:2, LogCarter1, double(1), 5, 2)
  lnTotal <- TreeTools::LnUnrooted(7)

  expect_equal(ms + lnTotal, lc, tolerance = 1e-12)
})

test_that("MaddisonSlatkin probabilities sum to 1 (2-state)", {
  for (a in 2:5) for (b in 2:a) {
    states <- c(as.integer(a), as.integer(b), 0L)
    maxSteps <- min(a, b)
    ms <- MaddisonSlatkin(1:maxSteps, states)
    expect_equal(sum(exp(ms)), 1, tolerance = 1e-10,
                 label = paste0("sum P for (", a, ",", b, ")"))
  }
})

test_that("MaddisonSlatkin probabilities sum to 1 (3-state)", {
  # (3, 2, 2) on 7 tips, min steps = 2
  states <- c(3L, 2L, 0L, 2L, 0L, 0L, 0L)
  n <- sum(states)
  ms <- MaddisonSlatkin(2:(n - 1L), states)
  expect_equal(sum(exp(ms[is.finite(ms)])), 1, tolerance = 1e-10)

  # (4, 3, 2) on 9 tips, min steps = 2
  states2 <- c(4L, 3L, 0L, 2L, 0L, 0L, 0L)
  n2 <- sum(states2)
  ms2 <- MaddisonSlatkin(2:(n2 - 1L), states2)
  expect_equal(sum(exp(ms2[is.finite(ms2)])), 1, tolerance = 1e-10)
})

test_that("MaddisonSlatkin probabilities sum to 1 (4-state)", {
  # (3, 2, 2, 2) on 9 tips, min steps = 3
  states <- integer(2^4 - 1)
  states[1] <- 3L  # state 1
  states[2] <- 2L  # state 2
  states[4] <- 2L  # state 3
  states[8] <- 2L  # state 4
  n <- sum(states)
  ms <- MaddisonSlatkin(3:(n - 1L), states)
  expect_equal(sum(exp(ms[is.finite(ms)])), 1, tolerance = 1e-10)
})

test_that("MaddisonSlatkin handles minimum step count correctly", {
  # 2-state: min steps = 1
  states <- c(4L, 3L, 0L)
  ms <- MaddisonSlatkin(1L, states)
  expect_true(is.finite(ms))
  expect_true(ms < 0)

  # 3-state: min steps = 2
  states3 <- c(3L, 3L, 0L, 3L, 0L, 0L, 0L)
  ms0 <- MaddisonSlatkin(1L, states3)
  expect_equal(ms0, -Inf)  # 1 step is impossible for 3 states

  ms2 <- MaddisonSlatkin(2L, states3)
  expect_true(is.finite(ms2))
  expect_true(ms2 < 0)
})

test_that("MaddisonSlatkin rejects invalid inputs", {
  expect_error(MaddisonSlatkin(1L, integer(0)))
  expect_error(MaddisonSlatkin(1L, c(-1L, 3L, 0L)))
})

test_that("MaddisonSlatkin_clear_cache runs without error", {
  expect_silent(MaddisonSlatkin_clear_cache())
})

test_that("MaddisonSlatkin with 5 states", {
  # (2,2,2,2,2) n=10: feasible (~0.6s); (3,2,2,2,2) n=11 blows up.
  states <- integer(31)
  states[1] <- 2L
  states[2] <- 2L
  states[4] <- 2L
  states[8] <- 2L
  states[16] <- 2L
  n <- sum(states)
  # min steps = 4 (one fewer than number of states)
  ms <- MaddisonSlatkin(4:(n - 1L), states)
  # On slow CI machines the 2s budget may be exceeded, yielding NA.
  # Skip the value checks in that case — the budget itself is correct
  # behaviour; we just can't verify the math when it fires.
  if (any(is.na(ms))) {
    skip("MaddisonSlatkin 5-state computation hit time budget")
  }
  expect_equal(sum(exp(ms[is.finite(ms)])), 1, tolerance = 1e-10)

  # Known value: (2,2,2,2,2) at 4 steps
  expect_equal(MaddisonSlatkin(4, states), -6.851185, tolerance = 1e-4)
})

test_that(".MSSplitCount is correct for known cases", {
  skip_on_cran()
  sc <- TreeSearch:::.MSSplitCount
  thresh <- TreeSearch:::.MS_SC_THRESHOLD

  # Balanced partitions at the empirical fast/blowup boundary
  expect_equal(sc(c(13L, 12L, 12L)), 133)         # k=3 n=37
  expect_equal(sc(c(13L, 13L, 12L)), 140)         # k=3 n=38
  expect_equal(sc(c(5L, 4L, 4L, 4L)),  95)        # k=4 n=17
  expect_equal(sc(c(5L, 5L, 4L, 4L)), 110)        # k=4 n=18
  expect_equal(sc(c(3L, 3L, 3L, 2L, 2L)),  96)    # k=5 n=13
  expect_equal(sc(c(3L, 3L, 3L, 3L, 2L)), 124)    # k=5 n=14

  # Skewed partitions are cheap even at large n
  expect_lt(sc(c(45L, 3L, 2L)), thresh[3])        # k=3 n=50
  expect_lt(sc(c(25L, 3L, 1L, 1L)), thresh[4])    # k=4 n=30
  expect_lt(sc(c(16L, 2L, 1L, 1L, 0L)), thresh[5]) # k=5 n=20

  # Base cases
  expect_equal(sc(integer(0)), 0)
  expect_equal(sc(c(1L, 1L)), 1)
  expect_equal(sc(c(0L, 3L)), 1)

  # Thresholds gate boundary cases correctly (correct bitmask encoding)
  expect_lte(sc(c(9L, 9L, 9L)),      thresh[3])   # k=3 n=27 sc=75: at limit
  expect_gt( sc(c(10L, 9L, 9L)),     thresh[3])   # k=3 n=28 sc=80: blocked
  expect_lte(sc(c(4L, 3L, 3L, 3L)),  thresh[4])   # k=4 n=13 sc=50: at limit
  expect_gt( sc(c(4L, 4L, 3L, 3L)),  thresh[4])   # k=4 n=14 sc=60: blocked
  expect_lte(sc(c(2L, 2L, 2L, 2L, 1L)), thresh[5]) # k=5 n=9  sc=35: at limit
  expect_gt( sc(c(2L, 2L, 2L, 2L, 2L)), thresh[5]) # k=5 n=10 sc=51: blocked
})


test_that("StepInformation() falls back instead of hanging when the exact memo cache overflows", {
  # A skewed 3-state character on many tips has a low split-count, so the
  # feasibility gate classes it feasible for the exact MaddisonSlatkin recursion
  # -- but on this many tips it generates more distinct (token, leaves) keys than
  # the solver's fixed-capacity memo tables reserve.  Before the capacity guard,
  # once such a table filled, probe_slot() spun forever (a hang that the 2 s
  # wall-clock budget could not interrupt, because the loop is outside
  # LogB/LogPVec).  On x86-64 the time budget usually won the race and bailed to
  # Monte Carlo; on Apple Silicon the table filled first and the search hung
  # (observed as a 6 h --run-donttest CI timeout).  The solver must now detect
  # the impending overflow and fall back to the MC approximation instead.
  #
  # Two guards can trigger that fallback, and which fires is a timing-dependent
  # race: the capacity guard once a memo table hits its reserved size, or the
  # 2 s wall-clock budget.  A fast machine reaches the capacity first (cache
  # warning); a slow one reaches 2 s first (time-budget warning).  Both warnings
  # contain "exceeded" and both are correct outcomes -- the only wrong outcome is
  # a hang.  So assert completion with finite values (the anti-hang property) and
  # that *some* fallback warning fired, without pinning which guard won the race.
  char <- rep(c("0", "1", "2"), c(42L, 9L, 2L))  # == inapplicable Agnarsson2004 col 83
  warnings_seen <- character(0)
  si <- withCallingHandlers(
    StepInformation(char, n_mc = 1000L),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_type(si, "double")
  expect_true(length(si) >= 1L && all(is.finite(si)))
  expect_match(paste(warnings_seen, collapse = "\n"), "exceeded")
})
