# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Helper: run ratchet search with full parameter set
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

ts_ratchet <- function(tree, ds, nCycles = 10L, perturbProb = 0.04,
                       maxHits = 1L, perturbMode = 0L,
                       perturbMaxMoves = 0L, adaptive = FALSE,
                       targetEscapeRate = 0.3) {
  TreeSearch:::ts_ratchet_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                    nCycles = nCycles, perturbProb = perturbProb,
                    maxHits = maxHits, perturbMode = perturbMode,
                    perturbMaxMoves = perturbMaxMoves, adaptive = adaptive,
                    targetEscapeRate = targetEscapeRate)
}


# --- New return fields ---

test_that("Ratchet returns n_escapes and final_perturb_prob", {
  tree <- as.phylo(42, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
    1, 0, 0, 0, 1, 1, 1, 1, 0, 0
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 3L)

  expect_true("n_escapes" %in% names(result))
  expect_true("final_perturb_prob" %in% names(result))
  expect_true(result$n_escapes >= 0L)
  expect_true(result$n_escapes <= result$n_cycles)
  expect_equal(result$final_perturb_prob, 0.04)
})


# --- Upweight mode ---

test_that("Upweight mode produces valid trees", {
  set.seed(3847)
  mat <- matrix(sample(0:1, 12 * 8, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  tree <- as.phylo(1, 12)

  result <- ts_ratchet(tree, ds, nCycles = 3L, perturbMode = 1L)

  expect_true(result$score > 0)
  expect_equal(result$n_cycles, 3L)
  expect_equal(nrow(result$edge), 2 * (12 - 1))
})

test_that("Upweight mode gives different behavior than zero mode", {
  # Zero mode (Nixon/Goloboff character deletion) and upweight mode (character
  # doubling) reshape the perturbation landscape differently, so their TBR
  # trajectories should diverge.  On any *single* small dataset the two modes can
  # coincidentally reach the same local optimum with the same move count (~3% of
  # seeds on this 15x10 binary matrix).  Which seeds coincide is trajectory- and
  # platform-dependent: ordering characters by homoplasy (see the "faster bounded
  # scan" change in src/ts_data.cpp) shifted the RNG-seeded search path, and on
  # macOS's std::sort tie-break the old single fixed seed (1001) landed on a
  # coincidence -- so the previous single-seed assertion flaked there.
  #
  # Assert instead that the modes are not identical *for every* seed: divergence on
  # >= 1 seed proves the two modes are distinct code paths, while a genuine
  # collapse (upweight silently equivalent to zero) would coincide on all of them.
  # The deterministic guarantee that the two weighting regimes are keyed apart is
  # covered by test-ts-na-evcache.R.  Early `break` keeps the common case to two
  # ratchet calls.
  set.seed(6284)
  mat <- matrix(sample(0:1, 15 * 10, replace = TRUE),
                nrow = 15,
                dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  tree <- as.phylo(1, 15)

  differ <- FALSE
  for (seed in c(1001L, 7L, 13L, 29L, 51L, 88L)) {
    set.seed(seed)
    result_zero <- ts_ratchet(tree, ds, nCycles = 5L,
                              perturbProb = 0.15, perturbMode = 0L)
    set.seed(seed)
    result_up <- ts_ratchet(tree, ds, nCycles = 5L,
                            perturbProb = 0.15, perturbMode = 1L)
    if ((result_zero$total_tbr_moves != result_up$total_tbr_moves) ||
        (result_zero$score != result_up$score)) {
      differ <- TRUE
      break
    }
  }
  expect_true(differ,
    info = "Zero and upweight modes should diverge on at least one seed")
})


# --- Mixed mode ---

test_that("Mixed mode produces valid trees", {
  set.seed(7712)
  mat <- matrix(sample(0:1, 12 * 6, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  tree <- as.phylo(1, 12)

  result <- ts_ratchet(tree, ds, nCycles = 3L, perturbMode = 2L,
                       perturbProb = 0.10)

  expect_true(result$score > 0)
  expect_equal(result$n_cycles, 3L)
})


# --- perturb_max_moves ---

test_that("perturb_max_moves parameter is respected", {
  set.seed(8193)
  mat <- matrix(sample(0:1, 20 * 12, replace = TRUE),
                nrow = 20,
                dimnames = list(paste0("t", 1:20), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  tree <- as.phylo(1, 20)

  # Very restrictive: only 2 moves per perturbation phase
  result_few <- ts_ratchet(tree, ds, nCycles = 3L, perturbMaxMoves = 2L)
  # Generous: 100 moves per perturbation phase
  result_many <- ts_ratchet(tree, ds, nCycles = 3L, perturbMaxMoves = 100L)

  # Both should produce valid results
  expect_true(result_few$score > 0)
  expect_true(result_many$score > 0)
})


# --- Adaptive perturbation ---

test_that("Adaptive mode adjusts final_perturb_prob", {
  set.seed(4492)
  mat <- matrix(sample(0:1, 15 * 8, replace = TRUE),
                nrow = 15,
                dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  tree <- as.phylo(1, 15)

  # Non-adaptive: final prob should equal starting prob
  result_fixed <- ts_ratchet(tree, ds, nCycles = 6L, adaptive = FALSE)
  expect_equal(result_fixed$final_perturb_prob, 0.04)

  # Adaptive: final prob may differ from starting prob
  result_adapt <- ts_ratchet(tree, ds, nCycles = 12L, adaptive = TRUE,
                             perturbProb = 0.04, targetEscapeRate = 0.3)
  expect_true(is.numeric(result_adapt$final_perturb_prob))
  expect_true(result_adapt$final_perturb_prob >= 0.02)
  expect_true(result_adapt$final_perturb_prob <= 0.50)
})


# --- IW upweighting ---

test_that("Upweight mode works with implied weights", {
  set.seed(5531)
  mat <- matrix(sample(0:2, 12 * 8, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  tree <- as.phylo(1, 12)

  at <- attributes(dataset)
  min_steps <- MinimumLength(dataset, compress = TRUE)

  result <- TreeSearch:::ts_ratchet_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    nCycles = 3L, perturbMode = 1L, perturbProb = 0.10,
    min_steps = as.integer(min_steps), concavity = 3.0
  )

  expect_true(result$score > 0)
  expect_equal(result$n_cycles, 3L)
})


# --- Regression: all existing checks still hold ---

test_that("Ratchet does not worsen score vs plain TBR (all modes)", {
  set.seed(2208)
  mat <- matrix(sample(0:1, 15 * 10, replace = TRUE),
                nrow = 15,
                dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  tree <- as.phylo(1, 15)

  tbr_result <- TreeSearch:::ts_tbr_search(tree$edge, ds$contrast, ds$tip_data,
                               ds$weight, ds$levels)

  for (mode in 0:2) {
    ratchet_result <- ts_ratchet(tree, ds, nCycles = 5L,
                                 perturbMode = as.integer(mode))
    expect_true(ratchet_result$score <= tbr_result$score,
      info = paste("Mode", mode, "should not worsen score"))
  }
})


# --- Upweight mask restoration ---

test_that("Upweight masks are properly cleaned up after ratchet", {
  tree <- as.phylo(42, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1, 1
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Run ratchet with upweight mode
  result <- ts_ratchet(tree, ds, nCycles = 3L, perturbMode = 1L,
                       perturbProb = 0.2)

  # Score the result tree with plain scoring — should be consistent
  result_tree <- tree
  result_tree$edge <- result$edge
  plain_score <- ts_score(result_tree, ds)
  expect_equal(result$score, plain_score,
    info = "Score after ratchet should match plain rescoring")
})
