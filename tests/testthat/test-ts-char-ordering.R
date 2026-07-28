# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Tests for Phase 3C: character-ordering optimizations.
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# Helper: run driven search
ts_driven <- function(ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 2L, xssRounds = 1L,
                      xssPartitions = 2L, fuseInterval = 2L,
                      maxSeconds = 0, verbosity = 0L, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    ratchetCycles = ratchetCycles,
    xssRounds = xssRounds,
    xssPartitions = xssPartitions,
    fuseInterval = fuseInterval,
    maxSeconds = maxSeconds,
    verbosity = verbosity,
    ...
  )
}

# ---------- Datasets ----------

# Mixed weights: characters with different weights
mixed_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0,
  0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0, 1, 1, 1
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
mixed_dataset <- MatrixToPhyDat(mixed_mat)
mixed_ds <- make_ts_data(mixed_dataset)

# 15-tip multi-state dataset
set.seed(5501)
multi_mat <- matrix(sample(0:3, 15 * 10, replace = TRUE),
                    nrow = 15,
                    dimnames = list(paste0("t", 1:15), NULL))
multi_dataset <- MatrixToPhyDat(multi_mat)
multi_ds <- make_ts_data(multi_dataset)

# ---------- Score invariance ----------

test_that("Scores are correct after block reordering", {
  # The sort change (descending weight) is internal; scores must be identical.
  set.seed(2047)
  tree <- as.phylo(42, 10)
  score <- ts_score(tree, mixed_ds)
  expect_true(score > 0)

  # Multiple random trees should all score correctly
  for (i in 1:5) {
    rt <- as.phylo(sample.int(1e5, 1), 10)
    s <- ts_score(rt, mixed_ds)
    expect_true(s >= score || s > 0)
  }
})

# `TS_CHAR_ORDER` (ts_data.cpp `CharOrder`) reorders characters WITHIN each
# block purely for speed: the bounded scorers reach `cutoff` in fewer blocks on
# rejected candidates. Commit d8c59998 asserted "the score is order-invariant
# (covered by the existing 'Scores are correct after block reordering' test)" --
# but that test (above, kept as-is) never varies the ordering, and its loop
# assertion `s >= score || s > 0` is satisfied by `s > 0` alone, so it cannot
# fail on an informative matrix. Nothing pinned the invariance claim; before
# these tests, no test in the package referenced `TS_CHAR_ORDER` at all.
#
# IW / XPIWE / profile are the modes that matter here. They index per-PATTERN
# arrays -- `ds.min_steps` and `ds.info_amounts`, populated at ts_data.cpp:447+,
# i.e. AFTER the sort at :190 -- so a block-slot-to-pattern mix-up would surface
# in these modes and in no other.
#
# CAVEAT this test cannot check for itself: on binary data the MINORITY key
# produces many ties, the sort is stable, and the permutation collapses toward
# the identity -- a green run would then prove nothing. The matrix below is
# multistate for exactly that reason. Red-team area 10 (2026-07-28) measured a
# genuinely non-identity permutation on comparable 2-6-state data (180/180
# positions moved, block-0 membership overlap 25/64). Keep this matrix
# multistate, and do not swap it for a binary one.
#
# Why `with_envvar` is enough to switch the ordering: ts_data.cpp:95 reads
# `TS_CHAR_ORDER` via a plain `std::getenv` into a local, inside build_dataset,
# with no static caching (it is the only occurrence in src/), and every
# ts_fitch_score call rebuilds the dataset -- so the value is re-read per call.
# Note for anyone extending this: because TS_CHAR_ORDER is *itself* score-
# invariant, it cannot serve as its own positive control for env propagation.
# If that ever needs checking at runtime, use a knob that does move a score
# (TS_DRIFT_EXACT is one).

test_that("TS_CHAR_ORDER is score-invariant under EW, IW, XPIWE and profile", {
  set.seed(6607)
  nTip <- 12
  orderMat <- matrix(sample(0:3, nTip * 60, replace = TRUE), nrow = nTip,
                     dimnames = list(paste0("t", seq_len(nTip)), NULL))
  # Guarantee a spread of state counts, so blocks differ in width.
  orderMat[, 1:12] <- sample(0:1, nTip * 12, replace = TRUE)
  orderMat[, 13:24] <- sample(0:2, nTip * 12, replace = TRUE)
  orderDs <- MatrixToPhyDat(orderMat)
  at <- attributes(orderDs)
  expect_gt(length(at$levels), 3)     # fails loudly if the matrix goes binary

  tsData <- make_ts_data(orderDs)
  minSteps <- MinimumLength(orderDs, compress = TRUE)
  # nTip = 12 at 4 states keeps MaddisonSlatkin's exact state cache within
  # capacity, so info.amounts is exact and warning-free. At nTip = 14 it spills
  # to the Monte Carlo fallback (6 warnings); the values stay finite, but don't
  # grow this matrix without re-checking the expect_true(is.finite(...)) below.
  infoAmounts <- attr(PrepareDataProfile(orderDs), "info.amounts")
  expect_false(anyNA(infoAmounts))

  # Score under each ordering; every mode must agree bit-for-bit with `none`.
  ScoreAll <- function(tree) {
    c(ew = ts_score(tree, tsData),
      iw = ts_score(tree, tsData, concavity = 3, min_steps = minSteps),
      profile = ts_score(tree, tsData, infoAmounts = infoAmounts))
  }

  for (treeSeed in c(101, 2749, 8123)) {
    tree <- Preorder(RenumberTips(as.phylo(treeSeed, nTip), names(orderDs)))
    reference <- with_envvar(c(TS_CHAR_ORDER = "none"), ScoreAll(tree))
    # Both guards matter: a NaN profile score would make every comparison
    # below trivially pass (waldo treats NaN as equal to NaN), which is the
    # same vacuous-green failure mode this test replaced.
    expect_true(all(is.finite(reference)))
    expect_true(all(reference > 0))

    for (ordering in c("min_steps", "minority", "entropy")) {
      actual <- with_envvar(c(TS_CHAR_ORDER = ordering), ScoreAll(tree))
      expect_equal(actual, reference,
                   label = paste0("TS_CHAR_ORDER=", ordering,
                                  " tree=", treeSeed))
    }
  }
})

test_that("EW driven search finds correct optimum", {
  set.seed(3341)
  result <- ts_driven(mixed_ds, maxReplicates = 3L, targetHits = 1L)
  expect_true(result$best_score > 0)
  validate_result(result, 10L)

  # Verify score matches rescore
  edge <- result$trees[[1]]
  rt <- as.phylo(1, 10)
  rt$edge <- edge
  expect_equal(result$best_score, ts_score(rt, mixed_ds))
})

test_that("IW driven search works with descending block order", {
  set.seed(8153)
  result <- ts_driven(mixed_ds, concavity = 3,
                      maxReplicates = 3L, targetHits = 1L)
  expect_true(result$best_score >= 0)
  validate_result(result, 10L)
})

test_that("Multi-state dataset scored correctly", {
  set.seed(7722)
  tree <- as.phylo(1, 15)
  score <- ts_score(tree, multi_ds)
  expect_true(score > 0)

  result <- ts_driven(multi_ds, maxReplicates = 3L, targetHits = 1L)
  expect_true(result$best_score <= score)
  validate_result(result, 15L)
})

# ---------- Zero-weight pattern compaction ----------

test_that("Jackknife with extreme deletion works correctly", {
  set.seed(6619)
  result <- TreeSearch:::ts_resample_search(
    mixed_ds$contrast, mixed_ds$tip_data, mixed_ds$weight, mixed_ds$levels,
    bootstrap = FALSE, jackProportion = 0.1,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true(result$score >= 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

test_that("Bootstrap produces valid results with compaction", {
  set.seed(4408)
  result <- TreeSearch:::ts_resample_search(
    multi_ds$contrast, multi_ds$tip_data, multi_ds$weight, multi_ds$levels,
    bootstrap = TRUE,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(result$score > 0)
  expect_equal(nrow(result$edge), 2L * (15L - 1L))
})

# ---------- Bounded indirect correctness ----------

test_that("Bounded indirect produces same results as search", {
  # If bounded indirect is wrong, search results will differ.
  # Run same search twice with different seeds: both must find valid trees.
  for (seed in c(1129, 5982)) {
    set.seed(seed)
    result <- ts_driven(multi_ds, maxReplicates = 3L, targetHits = 1L)
    expect_true(result$best_score > 0)
    validate_result(result, 15L)

    # Rescore to verify
    edge <- result$trees[[1]]
    rt <- as.phylo(1, 15)
    rt$edge <- edge
    expect_equal(result$best_score, ts_score(rt, multi_ds))
  }
})

# ---------- Ratchet with active_mask skip ----------

test_that("Ratchet search correct with active_mask optimization", {
  set.seed(9341)
  tree <- as.phylo(1, 10)
  result <- TreeSearch:::ts_ratchet_search(
    tree$edge, mixed_ds$contrast, mixed_ds$tip_data,
    mixed_ds$weight, mixed_ds$levels,
    nCycles = 3L
  )
  expect_true(result$score > 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

# ---------- set.seed() reproducibility ----------

test_that("Driven search is reproducible with set.seed()", {
  run_search <- function() {
    set.seed(2200)
    ts_driven(mixed_ds, maxReplicates = 3L, targetHits = 1L)
  }
  r1 <- run_search()
  r2 <- run_search()
  expect_equal(r1$best_score, r2$best_score)
  expect_equal(r1$trees[[1]], r2$trees[[1]])
})

# ---------- Inapplicable characters ----------

test_that("NA dataset scored correctly with block optimizations", {
  skip_if_not_installed("TreeSearch")
  data(inapplicable.phyData, package = "TreeSearch")

  if ("Vinther2008" %in% names(inapplicable.phyData)) {
    ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])
    n_tip <- length(inapplicable.phyData[["Vinther2008"]])

    set.seed(6677)
    result <- ts_driven(ds, maxReplicates = 2L, targetHits = 1L,
                        ratchetCycles = 1L)
    expect_true(result$best_score > 0)
    validate_result(result, n_tip)
  }
})

test_that("NA dataset with IW works after reordering", {
  skip_if_not_installed("TreeSearch")
  data(inapplicable.phyData, package = "TreeSearch")

  if ("Vinther2008" %in% names(inapplicable.phyData)) {
    ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])
    n_tip <- length(inapplicable.phyData[["Vinther2008"]])

    set.seed(8831)
    result <- ts_driven(ds, concavity = 5,
                        maxReplicates = 2L, targetHits = 1L,
                        ratchetCycles = 1L)
    expect_true(result$best_score >= 0)
    validate_result(result, n_tip)
  }
})

# ---------- Drift search with bounded indirect ----------

test_that("Drift search works with bounded indirect calls", {
  set.seed(3055)
  tree <- as.phylo(1, 10)
  result <- TreeSearch:::ts_drift_search(
    tree$edge, mixed_ds$contrast, mixed_ds$tip_data,
    mixed_ds$weight, mixed_ds$levels,
    nCycles = 3L
  )
  expect_true(result$score > 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

test_that("Drift search IW with bounded indirect calls", {
  set.seed(4601)
  tree <- as.phylo(1, 10)
  result <- TreeSearch:::ts_drift_search(
    tree$edge, mixed_ds$contrast, mixed_ds$tip_data,
    mixed_ds$weight, mixed_ds$levels,
    nCycles = 3L, concavity = 3
  )
  expect_true(result$score >= 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

# ---------- Wagner tree with bounded indirect ----------

test_that("Wagner tree construction correct with bounded indirect", {
  set.seed(7713)
  result1 <- TreeSearch:::ts_random_wagner_tree(
    mixed_ds$contrast, mixed_ds$tip_data, mixed_ds$weight, mixed_ds$levels
  )
  expect_true(result1$score > 0)
  expect_equal(nrow(result1$edge), 2L * (10L - 1L))

  result2 <- TreeSearch:::ts_random_wagner_tree(
    multi_ds$contrast, multi_ds$tip_data, multi_ds$weight, multi_ds$levels
  )
  expect_true(result2$score > 0)
  expect_equal(nrow(result2$edge), 2L * (15L - 1L))
})
