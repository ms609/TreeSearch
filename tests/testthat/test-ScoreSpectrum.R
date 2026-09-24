# Tier 1 (CRAN): fast pure-R tests for ScoreSpectrum().
# No C++ calls — no skip_on_cran() needed.

test_that("ScoreSpectrum handles edge cases gracefully", {
  # Empty input
  sp0 <- ScoreSpectrum(numeric(0))
  expect_equal(sp0$n_replicates, 0L)
  expect_equal(sp0$observed_levels, 0L)
  expect_true(is.na(sp0$coverage))

  # Single replicate
  sp1 <- ScoreSpectrum(42)
  expect_equal(sp1$n_replicates, 1L)
  expect_equal(sp1$observed_levels, 1L)
  expect_true(is.na(sp1$coverage))

  # Non-finite values are stripped (NA, Inf, -Inf all removed)
  sp_na <- ScoreSpectrum(c(10, NA, Inf, -Inf, 20, 10))
  expect_equal(sp_na$n_replicates, 3L)  # 3 finite values: 10, 20, 10
})

test_that("ScoreSpectrum computes correct Chao1 and coverage estimates", {
  # 10 replicates: 3 score-A (once each), 2 score-B (twice each), 1 score-C (five times)
  # score A seen 1x, score B seen 2x ... let's construct explicit
  # scores: 100, 102, 104, 101, 101, 103, 103, 103, 103, 103
  # abundance: 100->1, 101->2, 102->1, 103->5, 104->1
  # S_obs = 5, n = 10
  # f1 = 3 (100, 102, 104 each seen once), f2 = 1 (101 seen twice)
  # Chao1 = 5 + 3^2/(2*1) = 5 + 4.5 = 9.5
  # coverage = 1 - 3/10 = 0.7
  scores <- c(100, 101, 101, 102, 103, 103, 103, 103, 103, 104)
  sp <- ScoreSpectrum(scores, tol = 0)

  expect_equal(sp$n_replicates, 10L)
  expect_equal(sp$observed_levels, 5L)
  expect_equal(sp$estimated_levels, 9.5)
  expect_equal(sp$coverage, 0.7)
  expect_equal(sp$best_score, 100)
  expect_equal(sp$best_score_reps, 1L)
  # Frequency spectrum
  expect_equal(sp$f[["1"]], 3L)
  expect_equal(sp$f[["2"]], 1L)
  expect_equal(sp$f[["5"]], 1L)
  expect_equal(sum(sp$f), 5L)   # total = S_obs
})

test_that("ScoreSpectrum uses bias-corrected form when f2 = 0", {
  # All scores are singletons: f1 = n, f2 = 0
  # Chao1_bc = S_obs + f1*(f1-1)/2 = 5 + 5*4/2 = 15
  scores <- c(1, 2, 3, 4, 5)
  sp <- ScoreSpectrum(scores, tol = 0)
  expect_equal(sp$observed_levels, 5L)
  expect_equal(sp$estimated_levels, 5 + 5 * 4 / 2)
  # coverage = 1 - 5/5 = 0
  expect_equal(sp$coverage, 0)
})

test_that("ScoreSpectrum returns coverage = 1 when all scores are identical", {
  # Single score value, many replicates: f1 = 0, S_obs = 1
  scores <- rep(57, 20)
  sp <- ScoreSpectrum(scores, tol = 0)
  expect_equal(sp$observed_levels, 1L)
  expect_equal(sp$estimated_levels, 1)
  expect_equal(sp$coverage, 1)
  expect_equal(sp$unseen_fraction, 0)
  expect_equal(sp$best_score_reps, 20L)
})

test_that("ScoreSpectrum bins floating-point scores with tolerance", {
  # Two scores that differ by 1e-5 should be treated as equal at default tol
  scores <- c(100.00000, 100.00001, 100.00001, 101.00000)
  sp <- ScoreSpectrum(scores, tol = 1e-4)
  # 100.00000 and 100.00001 both round to 100.0000 at tol=1e-4 -> same bin
  expect_equal(sp$observed_levels, 2L)

  # With tol = 0 they should be distinct
  sp_exact <- ScoreSpectrum(scores, tol = 0)
  expect_equal(sp_exact$observed_levels, 3L)
})

test_that("ScoreSpectrum accepts raw numeric vector or multiPhylo", {
  scores <- c(10, 10, 10, 20, 30)
  sp_vec <- ScoreSpectrum(scores)
  expect_s3_class(sp_vec, "ScoreSpectrum")

  # Simulate a multiPhylo with replicate_scores attribute
  fake_trees <- structure(list(), class = "multiPhylo")
  attr(fake_trees, "replicate_scores") <- scores
  sp_mp <- ScoreSpectrum(fake_trees)
  expect_equal(sp_mp$n_replicates, sp_vec$n_replicates)
  expect_equal(sp_mp$coverage, sp_vec$coverage)
})

test_that("ScoreSpectrum errors informatively on bad input", {
  expect_error(ScoreSpectrum("not numeric"), "numeric vector")

  # multiPhylo without replicate_scores
  fake <- structure(list(), class = "multiPhylo")
  expect_error(ScoreSpectrum(fake), "replicate_scores")
})

test_that("print.ScoreSpectrum is callable without error", {
  sp <- ScoreSpectrum(c(100, 100, 101, 102, 100))
  expect_output(print(sp), "coverage")
  expect_output(print(sp), "Best score")
})

test_that("print.ScoreSpectrum handles insufficient replicates", {
  sp <- ScoreSpectrum(numeric(1))
  expect_output(print(sp), "insufficient")
})
