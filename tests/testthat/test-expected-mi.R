# Tier 1: arithmetic only, no search, whole file under a second.
#
# An independent reference for the expected mutual information under the
# hypergeometric null.  lchoose() works in log space, so unlike the C++
# recurrence it cannot underflow at the tails of the distribution.
ReferenceEmi <- function(ni, nj) {
  a <- ni[[1]]
  n <- sum(ni)
  emi <- 0
  for (mj in nj) {
    k <- max(0, a + mj - n):min(a, mj)
    logP <- lchoose(mj, k) + lchoose(n - mj, a - k) - lchoose(n, a)
    p <- exp(logP)
    kOut <- mj - k
    emi <- emi +
      sum(p * ifelse(k > 0, (k / n) * log2(k * n / (a * mj)), 0)) +
      sum(p * ifelse(kOut > 0, (kOut / n) * log2(kOut * n / ((n - a) * mj)), 0))
  }
  emi
}

test_that("expected_mi() is correct for large balanced partitions", {
  # P(K = kmin) is around 2^-1197 at N = 1200; a recurrence seeded there
  # returns exactly zero for every k.
  expect_equal(expected_mi(c(550L, 550L), c(550L, 550L)),
               ReferenceEmi(c(550L, 550L), c(550L, 550L)), tolerance = 1e-9)
  expect_equal(expected_mi(c(600L, 600L), c(600L, 600L)),
               ReferenceEmi(c(600L, 600L), c(600L, 600L)), tolerance = 1e-9)
  expect_equal(expected_mi(c(1000L, 1000L), c(1000L, 1000L)),
               ReferenceEmi(c(1000L, 1000L), c(1000L, 1000L)), tolerance = 1e-9)

  # `nj` as ClusteringConcordance() supplies it: a tabulate() over states,
  # which need not be two
  expect_equal(expected_mi(c(600L, 600L), c(300L, 300L, 300L, 300L)),
               ReferenceEmi(c(600L, 600L), c(300L, 300L, 300L, 300L)),
               tolerance = 1e-9)
  expect_equal(expected_mi(c(437L, 1063L), c(211L, 396L, 893L)),
               ReferenceEmi(c(437L, 1063L), c(211L, 396L, 893L)),
               tolerance = 1e-9)

  # The value the caller that motivated the fix actually receives
  expect_equal(TreeSearch:::.ExpectedMI(c(600L, 600L), c(600L, 600L)),
               ReferenceEmi(c(600L, 600L), c(600L, 600L)), tolerance = 1e-9)

  # Chance-corrected mutual information is positive and decreases with N
  balanced <- vapply(c(500L, 1000L, 1100L, 1200L, 2000L), function(n) {
    expected_mi(c(n %/% 2L, n %/% 2L), c(n %/% 2L, n %/% 2L))
  }, double(1))
  expect_true(all(balanced > 0))
  expect_true(all(diff(balanced) < 0))
})

test_that("expected_mi() is unchanged for small partitions", {
  # Values produced before the recurrence was re-anchored at the mode,
  # in the regime where seeding it at kmin was safe.
  expect_equal(expected_mi(c(3L, 4L), c(2L, 5L)),
               0.15383715015513183, tolerance = 1e-10)
  expect_equal(expected_mi(c(9L, 11L), c(4L, 7L, 9L)),
               0.084112593221791668, tolerance = 1e-10)
  expect_equal(expected_mi(c(50L, 50L), c(50L, 50L)),
               0.007323652940324857, tolerance = 1e-10)
  expect_equal(expected_mi(c(37L, 163L), c(11L, 60L, 129L)),
               0.0077921375666311615, tolerance = 1e-10)
  expect_equal(expected_mi(c(500L, 500L), c(500L, 500L)),
               0.00072243147032421289, tolerance = 1e-10)
  expect_equal(expected_mi(c(400L, 800L), c(300L, 400L, 500L)),
               0.0012047105794341356, tolerance = 1e-10)
  expect_equal(expected_mi(c(1L, 6L), c(3L, 4L)), 0.15809905413668374,
               tolerance = 1e-10)
  expect_equal(expected_mi(c(0L, 7L), c(3L, 4L)), 0)
  expect_equal(expected_mi(c(7L, 0L), c(3L, 4L)), 0)
})

test_that("expected_mi() rejects an `ni` that is not a pair", {
  expect_error(expected_mi(3L, c(2L, 5L)), "length 2")
  expect_error(expected_mi(integer(0), c(2L, 5L)), "length 2")
  expect_error(expected_mi(c(1L, 2L, 4L), c(2L, 5L)), "length 2")
})

test_that("expected_mi() agrees across the factorial lookup boundary", {
  # N exceeds the 8192-entry log-factorial table, so l2factorial() must
  # return matching values from the table and from its lgamma() fallback.
  expect_equal(expected_mi(c(4500L, 4500L), c(4500L, 4500L)),
               ReferenceEmi(c(4500L, 4500L), c(4500L, 4500L)),
               tolerance = 1e-8)
  expect_equal(expected_mi(c(3000L, 7000L), c(4096L, 5904L)),
               ReferenceEmi(c(3000L, 7000L), c(4096L, 5904L)),
               tolerance = 1e-8)
})

test_that("mi_key() distinguishes block sizes above 65535", {
  # Sorting is only sound because expected_mi() is invariant under both
  # canonicalizations the key applies.
  expect_equal(expected_mi(c(3L, 61L), c(30L, 31L)),
               expected_mi(c(61L, 3L), c(31L, 30L)))
  expect_identical(TreeSearch:::mi_key(c(3L, 61L), c(30L, 31L)),
                   TreeSearch:::mi_key(c(61L, 3L), c(31L, 30L)))

  # Fixed-width hex of the sorted values.  0x12345678 has eight distinct
  # nibbles, so an emit that mis-shifts one is caught here; injectivity
  # alone would not notice.
  expect_identical(
    TreeSearch:::mi_key(c(65597L, 3L), c(305419896L, 30L)),
    paste(sprintf("%08x", c(3L, 65597L, 30L, 305419896L)), collapse = "")
  )

  # Block sizes differing by a multiple of 65536 must not share a key
  aliases <- c(60L, 61L, 65596L, 65597L, 131133L)
  keys <- vapply(aliases, function(n) {
    TreeSearch:::mi_key(c(3L, n), c(30L, 31L))
  }, character(1))
  expect_equal(anyDuplicated(keys), 0L)

  # The cached value must belong to the partition asked for.  Populate the
  # small key first, so a colliding key would return it.
  expect_equal(TreeSearch:::.ExpectedMI(c(3L, 61L), c(30L, 31L)),
               expected_mi(c(3L, 61L), c(30L, 31L)))
  expect_equal(TreeSearch:::.ExpectedMI(c(3L, 65597L), c(30L, 31L)),
               expected_mi(c(3L, 65597L), c(30L, 31L)))
})

test_that("quartet_concordance() rejects negative state codes", {
  splits <- matrix(c(TRUE, TRUE, FALSE, FALSE), ncol = 1)
  characters <- matrix(c(1L, 1L, 2L, 2L), ncol = 1)
  counts <- TreeSearch:::quartet_concordance(splits, characters)
  expect_equal(counts[["concordant"]], matrix(1))
  expect_equal(counts[["decisive"]], matrix(1))

  negative <- matrix(c(1L, -1L, 2L, 2L), ncol = 1)
  expect_error(TreeSearch:::quartet_concordance(splits, negative),
               "non-negative")

  # NA marks the absence of a state, and is not a negative code: a taxon
  # scored NA counts as if it were not in the matrix at all.  Six taxa give
  # a quartet count that a dropped taxon could change, unlike three.
  sixSplits <- matrix(c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), ncol = 1)
  sixChars <- matrix(c(1L, 1L, 2L, 2L, 1L, 2L), ncol = 1)
  sixCounts <- TreeSearch:::quartet_concordance(sixSplits, sixChars)
  expect_equal(sixCounts[["concordant"]], matrix(1))
  expect_equal(sixCounts[["decisive"]], matrix(5))
  expect_equal(
    TreeSearch:::quartet_concordance(rbind(sixSplits, TRUE),
                                     rbind(sixChars, NA_integer_)),
    sixCounts)
})
