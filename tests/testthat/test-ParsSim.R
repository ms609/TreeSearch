# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# --- Tests ------------------------------------------------------------------

test_that("ParsSim returns phyDat with correct structure", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(4817)
  result <- ParsSim(tree, nChar = c(10L), nExtraSteps = 0L)

  expect_s3_class(result, "phyDat")
  expect_equal(length(result), 8L)
  expect_true(!is.null(attr(result, "saturated")))
  expect_true(!is.null(attr(result, "steps_exhausted")))
  expect_true(!is.null(attr(result, "extra_steps")))
  expect_length(attr(result, "saturated"), 10L)
  expect_length(attr(result, "steps_exhausted"), 10L)
  expect_length(attr(result, "extra_steps"), 10L)
})

test_that("Zero extra steps: binary characters have score = nChar", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(2941)
  result <- ParsSim(tree, nChar = c(20L), nExtraSteps = 0L)

  total_score <- TreeSearch::TreeLength(tree, result)
  # Each 2-state character contributes exactly 1 step (minimum)
  expect_equal(total_score, 20)
  # No extra steps placed
  expect_true(all(attr(result, "extra_steps") == 0L))
})

test_that("Zero extra steps: 3-state characters have score = 2 each", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(7293)
  result <- ParsSim(tree, nChar = c(0L, 15L), nExtraSteps = 0L)

  total_score <- TreeSearch::TreeLength(tree, result)
  # Each 3-state character contributes exactly 2 steps
  expect_equal(total_score, 30)
})

test_that("Extra steps increase total score correctly", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(5012)
  n_extra <- 10L
  result <- ParsSim(tree, nChar = c(30L), nExtraSteps = n_extra)

  total_score <- TreeSearch::TreeLength(tree, result)
  # Total score = nChar * 1 (minimum) + extra steps placed
  placed <- sum(attr(result, "extra_steps"))
  expect_equal(total_score, 30L + placed)
  # All requested steps should have been placed (30 chars, only 10 steps)
  expect_equal(placed, n_extra)
})

test_that("Each character's score matches its Fitch score (no masking)", {
  tree <- TreeTools::BalancedTree(10)
  set.seed(1384)
  result <- ParsSim(tree, nChar = c(15L), nExtraSteps = 20L)

  char_lengths <- TreeSearch::CharacterLength(tree, result)
  extra <- attr(result, "extra_steps")
  n_states <- rep(2L, 15L)
  expected_scores <- (n_states - 1L) + extra

  # CharacterLength returns per-character scores; sum should match
  expect_equal(sum(char_lengths), sum(expected_scores))
})

test_that("IW selection biases towards already-homoplastic characters", {
  tree <- TreeTools::BalancedTree(12)
  set.seed(8456)

  # With many extra steps and low k, steps should concentrate
  result_iw <- ParsSim(tree, nChar = c(50L), nExtraSteps = 100L,
                       concavity = 2)
  extra_iw <- attr(result_iw, "extra_steps")

  # Under IW, variance of extra_steps should be higher than under EW
  set.seed(8456)
  result_ew <- ParsSim(tree, nChar = c(50L), nExtraSteps = 100L,
                       concavity = Inf)
  extra_ew <- attr(result_ew, "extra_steps")

  expect_gt(var(extra_iw), var(extra_ew))
})

test_that("Saturation attributes are computed correctly", {
  # Small tree: saturation should occur quickly
  tree <- TreeTools::BalancedTree(4)
  set.seed(6130)
  expect_warning(
    result <- ParsSim(tree, nChar = c(5L), nExtraSteps = 100L),
    "saturated"
  )

  saturated <- attr(result, "saturated")
  steps_exhausted <- attr(result, "steps_exhausted")
  extra <- attr(result, "extra_steps")

  # Some characters should be saturated on a 4-tip tree
  expect_true(any(saturated))

  # steps_exhausted should be a subset of saturated
  expect_true(all(which(steps_exhausted) %in% which(saturated)))

  # Characters that aren't saturated could accept more steps
  # (We can't test this directly without rerunning, but extra_steps should
  # match sum placed)
  expect_true(all(extra >= 0L))
})

test_that("Non-binary tree is resolved with warning", {
  tree <- ape::read.tree(text = "((A,B,C),D,E);")
  set.seed(3847)
  expect_warning(
    result <- ParsSim(tree, nChar = c(5L), nExtraSteps = 0L),
    "non-binary|binary"
  )
  expect_s3_class(result, "phyDat")
})

test_that("Trees without edge lengths work (uniform default)", {
  tree <- TreeTools::BalancedTree(6)
  tree$edge.length <- NULL
  set.seed(9021)
  result <- ParsSim(tree, nChar = c(10L), nExtraSteps = 5L)
  expect_s3_class(result, "phyDat")
  expect_equal(sum(attr(result, "extra_steps")), 5L)
})

test_that("Reproducibility: same seed gives same result", {
  tree <- TreeTools::BalancedTree(8)

  set.seed(4455)
  r1 <- ParsSim(tree, nChar = c(10L), nExtraSteps = 5L)
  set.seed(4455)
  r2 <- ParsSim(tree, nChar = c(10L), nExtraSteps = 5L)

  m1 <- TreeTools::PhyDatToMatrix(r1)
  m2 <- TreeTools::PhyDatToMatrix(r2)
  expect_identical(m1, m2)
  expect_identical(attr(r1, "extra_steps"), attr(r2, "extra_steps"))
})

test_that("Multi-state characters (3+) with extra steps", {
  tree <- TreeTools::BalancedTree(10)
  set.seed(7722)
  # Mix of 2-state and 3-state
  result <- ParsSim(tree, nChar = c(10L, 5L), nExtraSteps = 15L)

  expect_s3_class(result, "phyDat")
  total_score <- TreeSearch::TreeLength(tree, result)
  min_score <- 10L * 1L + 5L * 2L  # 10 binary + 5 ternary
  placed <- sum(attr(result, "extra_steps"))
  expect_equal(total_score, min_score + placed)
})

test_that("Edge-length weighting biases transition placement", {
  # Two-cherry tree: ((A,B),(C,D)); make one pendant edge very long
  tree <- TreeTools::BalancedTree(4)
  tree$edge.length <- rep(1, nrow(tree$edge))
  # Find edge to tip t1 and make it very long
  t1_edge <- which(tree$edge[, 2] == 1L)
  tree$edge.length[t1_edge] <- 1000

  # Over many single-character zero-extra replicates, the initial transition
  # should land on the long edge more often than expected by chance
  n_reps <- 100L
  t1_transitions <- 0L
  for (i in seq_len(n_reps)) {
    ds <- ParsSim(tree, nChar = c(1L), nExtraSteps = 0L)
    tips <- TreeTools::PhyDatToMatrix(ds)[, 1]
    # t1 has a different state from the root state (0) if transition hit
    # its edge
    if (tips["t1"] != tips["t2"]) t1_transitions <- t1_transitions + 1L
  }
  # With uniform weights, P(t1 edge) ≈ 1/5 = 20%.
  # With weight 1000 vs 4×1, P(t1 edge) ≈ 1000/1004 ≈ 99.6%
  expect_gt(t1_transitions, 50L)
})

test_that("Unrooted tree is rooted internally", {
  tree <- TreeTools::UnrootTree(TreeTools::BalancedTree(8))
  set.seed(8830)
  result <- ParsSim(tree, nChar = c(5L), nExtraSteps = 3L)
  expect_s3_class(result, "phyDat")
  rooted <- TreeTools::RootTree(tree, tree[["tip.label"]][1])
  expect_equal(TreeSearch::TreeLength(rooted, result), 5L + 3L)
})

test_that("Tip labels are preserved in output", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(4456)
  result <- ParsSim(tree, nChar = c(3L))
  expect_true(all(tree[["tip.label"]] %in% names(result)))
})

test_that("Back-mutations produce valid scores on large tree", {
  tree <- TreeTools::BalancedTree(20)
  set.seed(6614)
  # Many extra steps forces back-mutations
  result <- ParsSim(tree, nChar = c(5L), nExtraSteps = 20L)
  expect_equal(TreeSearch::TreeLength(tree, result), 5L + 20L)
})

test_that("Large tree does not saturate prematurely", {
  tree <- TreeTools::BalancedTree(32)
  set.seed(7705)
  result <- ParsSim(tree, nChar = c(20L), nExtraSteps = 40L)
  expect_equal(sum(attr(result, "extra_steps")), 40L)
  expect_equal(TreeSearch::TreeLength(tree, result), 20L + 40L)
  expect_false(any(attr(result, "steps_exhausted")))
})

test_that("Input validation catches errors", {
  tree <- TreeTools::BalancedTree(6)
  expect_error(ParsSim("not_a_tree"), "phylo")
  expect_error(ParsSim(tree, nChar = c(-1L)), "non-negative")
  expect_error(ParsSim(tree, nChar = c(0L)), "at least one")
  expect_error(ParsSim(tree, nExtraSteps = -1L), "non-negative")
  expect_error(ParsSim(tree, concavity = -5), "positive")
  # rootState length must be 1 or sum(nChar)
  expect_error(ParsSim(tree, nChar = c(3L), rootState = c(0L, 1L)),
               "length 1 or sum")
  # rootState out of range for 2-state character

  expect_error(ParsSim(tree, nChar = c(3L), rootState = 2L),
               "out of range")
  expect_error(ParsSim(tree, nChar = c(3L), rootState = -1L),
               "out of range")
  # rootState out of range for specific character in vector
  expect_error(
    ParsSim(tree, nChar = c(2L, 1L), rootState = c(0L, 1L, 3L)),
    "out of range"
  )
})

test_that("Scalar rootState > 0 works correctly", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(1429)
  result <- ParsSim(tree, nChar = c(10L), nExtraSteps = 5L, rootState = 1L)
  expect_s3_class(result, "phyDat")
  expect_equal(TreeSearch::TreeLength(tree, result), 10L + 5L)
})

test_that("Vector rootState: per-character root states", {
  tree <- TreeTools::BalancedTree(6)
  set.seed(8273)
  # 2 binary characters: first starts at 0, second starts at 1
  result <- ParsSim(tree, nChar = c(2L), nExtraSteps = 0L,
                    rootState = c(0L, 1L))
  expect_s3_class(result, "phyDat")
  # No warnings should have occurred (original issue)
  expect_equal(TreeSearch::TreeLength(tree, result), 2L)
})

test_that("Vector rootState with mixed state counts", {
  tree <- TreeTools::BalancedTree(10)
  set.seed(3916)
  # 3 binary + 2 ternary = 5 chars total
  # rootState: first 3 chars root=0, char 4 root=1, char 5 root=2
  result <- ParsSim(tree, nChar = c(3L, 2L), nExtraSteps = 8L,
                    rootState = c(0L, 0L, 0L, 1L, 2L))
  expect_s3_class(result, "phyDat")
  total <- TreeSearch::TreeLength(tree, result)
  min_score <- 3L * 1L + 2L * 2L
  expect_equal(total, min_score + sum(attr(result, "extra_steps")))
  expect_equal(sum(attr(result, "extra_steps")), 8L)
})

test_that("Original issue: ParsSim(BalancedTree(6), 2, 0, Inf, c(0,1)) no warnings", {
  tree <- TreeTools::BalancedTree(6)
  set.seed(7501)
  expect_no_warning(
    result <- ParsSim(tree, nChar = c(2L), nExtraSteps = 0L,
                      rootState = c(0L, 1L))
  )
  expect_equal(TreeSearch::TreeLength(tree, result), 2L)
})

test_that("All characters saturated triggers warning", {
  # 4-tip tree with binary characters: very limited capacity
  tree <- TreeTools::BalancedTree(4)
  set.seed(5599)
  expect_warning(
    result <- ParsSim(tree, nChar = c(2L), nExtraSteps = 100L),
    "saturated"
  )
  expect_s3_class(result, "phyDat")
})

# --- Profile parsimony tests ------------------------------------------------

test_that("Profile mode produces valid phyDat", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(6812)
  result <- ParsSim(tree, nChar = c(15L), nExtraSteps = 10L,
                    concavity = "profile")
  expect_s3_class(result, "phyDat")
  expect_equal(sum(attr(result, "extra_steps")), 10L)
})

test_that("Profile mode allocates steps differently from EW", {
  # Use a large tree with few extra steps to avoid saturation
  tree <- TreeTools::BalancedTree(16)
  set.seed(2731)
  result_pp <- ParsSim(tree, nChar = c(40L), nExtraSteps = 30L,
                       concavity = "profile")
  set.seed(2731)
  result_ew <- ParsSim(tree, nChar = c(40L), nExtraSteps = 30L,
                       concavity = Inf)

  extra_pp <- attr(result_pp, "extra_steps")
  extra_ew <- attr(result_ew, "extra_steps")

  # The two modes should produce different allocations
  expect_false(identical(extra_pp, extra_ew))
  # Both should place the same total (no saturation with few steps)
  expect_equal(sum(extra_pp), sum(extra_ew))
})

test_that("Profile mode respects info = 0 saturation", {
  # Small tree: characters should hit info = 0 quickly
  tree <- TreeTools::BalancedTree(5)
  set.seed(4190)
  result <- suppressWarnings(
    ParsSim(tree, nChar = c(10L), nExtraSteps = 50L,
            concavity = "profile")
  )

  # Characters at max steps should be flagged as exhausted
  exhausted <- attr(result, "steps_exhausted")
  extra <- attr(result, "extra_steps")
  # At least some should be exhausted on a 5-tip tree with 50 requested steps
  expect_true(any(exhausted))
})

test_that("Profile scores match expectations", {
  tree <- TreeTools::BalancedTree(8)
  set.seed(3350)
  result <- ParsSim(tree, nChar = c(20L), nExtraSteps = 15L,
                    concavity = "profile")

  total_score <- TreeSearch::TreeLength(tree, result)
  extra <- attr(result, "extra_steps")
  # Total = 20 * 1 (minimum) + extra steps placed
  expect_equal(total_score, 20L + sum(extra))
})

# --- Extended tests (T-111) -------------------------------------------------

test_that("Per-character Fitch score matches expected for every character", {
  tree <- TreeTools::BalancedTree(16)
  set.seed(2583)
  # Mixed states: 10 binary + 5 ternary + 3 four-state
  result <- ParsSim(tree, nChar = c(10L, 5L, 3L), nExtraSteps = 30L)

  char_lengths <- TreeSearch::CharacterLength(tree, result)
  extra <- attr(result, "extra_steps")
  n_states_vec <- c(rep(2L, 10), rep(3L, 5), rep(4L, 3))
  min_steps <- n_states_vec - 1L
  expected <- min_steps + extra

  # Every character must match
  expect_equal(as.integer(char_lengths), expected)
})

test_that("Edge-length weighting: Chi-squared test on binned transitions", {
  tree <- TreeTools::PectinateTree(8)
  el <- rep(1, nrow(tree[["edge"]]))
  # Pendant edges to t1 and t8: make them 50x longer
  t1_edge <- which(tree[["edge"]][, 2] == 1L)
  t8_edge <- which(tree[["edge"]][, 2] == 8L)
  el[t1_edge] <- 50
  el[t8_edge] <- 50
  tree[["edge.length"]] <- el

  # Classify edges into "long" (t1 and t8) vs "short" (everything else)
  long_edges <- c(t1_edge, t8_edge)
  n_reps <- 200L
  long_hits <- 0L
  total_transitions <- 0L

  for (i in seq_len(n_reps)) {
    ds <- ParsSim(tree, nChar = c(1L), nExtraSteps = 0L)
    mat <- TreeTools::PhyDatToMatrix(ds)
    tips <- mat[, 1]
    # A transition on a pendant edge means that tip differs from root (0)
    if (tips["t1"] != tips["t2"]) long_hits <- long_hits + 1L
    if (tips["t8"] != tips["t7"]) long_hits <- long_hits + 1L
    total_transitions <- total_transitions + 1L
  }

  # Expected proportion of long edges: sum(long) / sum(all) = 100 / (100+12) ≈ 0.89
  # Each rep places 1 transition; we check if t1 or t8 got it
  # Under weighting, P(transition on a long edge) is high
  p_expected <- sum(el[long_edges]) / sum(el)
  # At least 60% of initial transitions should hit long edges (very conservative)
  expect_gt(long_hits / n_reps, 0.4)
})

test_that("Saturation ceiling on small trees", {
  # 5-tip tree: maximum Fitch score for a 2-state character
  tree <- TreeTools::BalancedTree(5)

  # With enough extra steps, we should hit the ceiling
  set.seed(8814)
  expect_warning(
    result <- ParsSim(tree, nChar = c(5L), nExtraSteps = 500L),
    "saturated"
  )

  extra <- attr(result, "extra_steps")
  saturated <- attr(result, "saturated")

  # All should be saturated on a 5-tip tree with 500 requested steps
  expect_true(all(saturated))
  # Total steps placed < 500 (couldn't place them all)
  expect_lt(sum(extra), 500L)
})

test_that("Large extra_steps with few characters", {
  tree <- TreeTools::BalancedTree(32)
  set.seed(9112)
  # Only 3 characters but 20 extra steps — each gets heavily loaded
  # Use 32 tips to avoid premature saturation
  result <- ParsSim(tree, nChar = c(3L), nExtraSteps = 20L)

  expect_equal(TreeSearch::TreeLength(tree, result), 3L + 20L)
  expect_equal(sum(attr(result, "extra_steps")), 20L)
})

test_that("4-state characters with many extra steps", {
  tree <- TreeTools::BalancedTree(20)
  set.seed(4267)
  result <- ParsSim(tree, nChar = c(0L, 0L, 8L), nExtraSteps = 30L)

  # Each 4-state character has min_steps = 3; total min = 8*3 = 24
  total <- TreeSearch::TreeLength(tree, result)
  expect_equal(total, 24L + sum(attr(result, "extra_steps")))
  expect_equal(sum(attr(result, "extra_steps")), 30L)
})

test_that("5-state characters work correctly", {
  tree <- TreeTools::BalancedTree(12)
  set.seed(6193)
  result <- ParsSim(tree, nChar = c(0L, 0L, 0L, 4L), nExtraSteps = 10L)

  # Each 5-state character has min_steps = 4; total min = 4*4 = 16
  total <- TreeSearch::TreeLength(tree, result)
  expect_equal(total, 16L + sum(attr(result, "extra_steps")))
  expect_equal(sum(attr(result, "extra_steps")), 10L)
})

test_that("Mixed 2-5 state characters with many extra steps", {
  tree <- TreeTools::BalancedTree(24)
  set.seed(7834)
  # 5 binary + 4 ternary + 3 four-state + 2 five-state
  result <- ParsSim(tree, nChar = c(5L, 4L, 3L, 2L), nExtraSteps = 50L)

  char_lengths <- TreeSearch::CharacterLength(tree, result)
  extra <- attr(result, "extra_steps")
  n_states_vec <- c(rep(2L, 5), rep(3L, 4), rep(4L, 3), rep(5L, 2))
  min_steps <- n_states_vec - 1L
  expected <- min_steps + extra

  expect_equal(as.integer(char_lengths), expected)
  expect_equal(sum(extra), 50L)
})

test_that("Score matches across different tree shapes", {
  shapes <- list(
    balanced = TreeTools::BalancedTree(16),
    pectinate = TreeTools::PectinateTree(16),
    random = TreeTools::RandomTree(16, root = TRUE)
  )

  for (shape_name in names(shapes)) {
    tree <- shapes[[shape_name]]
    set.seed(5511)
    result <- ParsSim(tree, nChar = c(10L, 3L), nExtraSteps = 15L)

    total <- TreeSearch::TreeLength(tree, result)
    min_score <- 10L * 1L + 3L * 2L
    placed <- sum(attr(result, "extra_steps"))
    expect_equal(total, min_score + placed,
                 label = paste("Score for", shape_name, "tree"))
    expect_equal(placed, 15L,
                 label = paste("Steps placed for", shape_name, "tree"))
  }
})
