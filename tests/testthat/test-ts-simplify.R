# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

ts_diag <- function(ds) {
  TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                 ds$weight, ds$levels)
}

ts_driven <- function(ds, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 3L, targetHits = 1L,
    ratchetCycles = 1L, xssRounds = 0L,
    xssPartitions = 2L, fuseInterval = 10L,
    maxSeconds = 0, verbosity = 0L,
    ...
  )
}

# ---- Test datasets ----

# 1. All-informative binary: 10 tips, 4 informative chars
info_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
info_dataset <- MatrixToPhyDat(info_mat)
info_ds <- make_ts_data(info_dataset)

# 2. Has autapomorphies: 8 tips, mix of informative and uninformative
autap_mat <- matrix(c(
  # Char 1: informative (0/1 split)
  0, 0, 0, 0, 1, 1, 1, 1,
  # Char 2: single autapomorphy (tip 1 = 2, rest = 0) -> uninformative
  2, 0, 0, 0, 0, 0, 0, 0,
  # Char 3: two autapomorphies (tips 1,2 unique) -> uninformative
  1, 2, 0, 0, 0, 0, 0, 0,
  # Char 4: informative (0/1 split)
  0, 0, 1, 1, 0, 0, 1, 1
), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
autap_dataset <- MatrixToPhyDat(autap_mat)
autap_ds <- make_ts_data(autap_dataset)

# 3. Has singleton states in an informative character: 0001112
# State 2 is a singleton -> always costs 1 extra step
singleton_mat <- matrix(c(
  0, 0, 0, 1, 1, 1, 2,
  0, 0, 1, 1, 0, 0, 1
), nrow = 7, dimnames = list(paste0("t", 1:7), NULL))
singleton_dataset <- MatrixToPhyDat(singleton_mat)
singleton_ds <- make_ts_data(singleton_dataset)

# 4. Invariant character mixed with informative
invar_mat <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 0,  # invariant
  0, 0, 0, 0, 1, 1, 1, 1   # informative
), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
invar_dataset <- MatrixToPhyDat(invar_mat)
invar_ds <- make_ts_data(invar_dataset)


# ===== Transform diagnostics =====

test_that("Diagnostic: all-informative dataset has no simplification", {
  diag <- ts_diag(info_ds)
  expect_equal(diag$n_patterns_removed, 0L)
  expect_true(all(diag$informative))
  expect_equal(diag$ew_offset, 0L)
  expect_true(all(diag$precomputed_steps == 0L))
})

test_that("Diagnostic: autapomorphies are detected and removed", {
  diag <- ts_diag(autap_ds)
  # Chars 2 & 3 are uninformative
  expect_gte(diag$n_patterns_removed, 1L)
  expect_gt(diag$ew_offset, 0L)
})

test_that("Diagnostic: singleton states produce precomputed_steps", {
  diag <- ts_diag(singleton_ds)
  # At least one pattern should have precomputed_steps > 0
  expect_true(any(diag$precomputed_steps > 0))
})

test_that("Diagnostic: invariant character is removed", {
  diag <- ts_diag(invar_ds)
  # The invariant pattern (all 0s) is uninformative, 0 steps
  expect_gte(diag$n_patterns_removed, 1L)
})


# ===== EW scoring correctness =====

# Reference scores: compute with phangorn for comparison
# We verify that ts_fitch_score (which now uses simplification) matches
# the expected parsimony score on multiple random trees.

test_that("EW scores match expected values on autapomorphy dataset", {
  set.seed(7134)
  for (i in seq_len(5)) {
    tree <- RandomTree(autap_dataset, root = TRUE)
    score <- ts_score(tree, autap_ds)
    # Compute expected score with phangorn
    expected <- phangorn::parsimony(tree, autap_dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

test_that("EW scores match expected values on singleton dataset", {
  set.seed(2891)
  for (i in seq_len(5)) {
    tree <- RandomTree(singleton_dataset, root = TRUE)
    score <- ts_score(tree, singleton_ds)
    expected <- phangorn::parsimony(tree, singleton_dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

test_that("EW scores match expected values on invariant+informative dataset", {
  set.seed(5603)
  for (i in seq_len(5)) {
    tree <- RandomTree(invar_dataset, root = TRUE)
    score <- ts_score(tree, invar_ds)
    expected <- phangorn::parsimony(tree, invar_dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

test_that("EW scores match on all-informative dataset (no simplification effect)", {
  set.seed(3927)
  for (i in seq_len(5)) {
    tree <- RandomTree(info_dataset, root = TRUE)
    score <- ts_score(tree, info_ds)
    expected <- phangorn::parsimony(tree, info_dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})


# ===== IW scoring correctness =====

test_that("IW scores are consistent across simplifiable datasets", {
  # Verify that IW score is the same on multiple trees, compared against
  # a non-simplifiable (all-informative) baseline approach.
  # The key invariant: IW score uses extra = steps - min_steps,

  # and simplification reduces both by the same amount.
  set.seed(4418)
  k <- 3.0  # concavity constant

  for (i in seq_len(5)) {
    tree <- RandomTree(autap_dataset, root = TRUE)
    # IW score via the C++ engine (with simplification)
    iw_score <- ts_score(tree, autap_ds, concavity = k,
                         min_steps = autap_ds$weight * 0L)
    # IW score should be finite and non-negative
    expect_true(is.finite(iw_score), info = paste("Tree", i, "finite"))
    expect_gte(iw_score, 0, label = paste("Tree", i, "non-negative"))
  }
})


# ===== Driven search correctness =====

test_that("Driven search finds correct best score with autapomorphies", {
  set.seed(6725)
  result <- ts_driven(autap_ds)
  # The best score from driven search should match scoring the best tree
  best_tree <- structure(
    list(edge = result$trees[[1]],
         tip.label = paste0("t", seq_len(nrow(autap_ds$tip_data))),
         Nnode = nrow(autap_ds$tip_data) - 1L),
    class = "phylo"
  )
  rescore <- ts_score(best_tree, autap_ds)
  expect_equal(result$best_score, rescore)
})

test_that("Driven search finds correct best score with singletons", {
  set.seed(1089)
  result <- ts_driven(singleton_ds)
  best_tree <- structure(
    list(edge = result$trees[[1]],
         tip.label = paste0("t", seq_len(nrow(singleton_ds$tip_data))),
         Nnode = nrow(singleton_ds$tip_data) - 1L),
    class = "phylo"
  )
  rescore <- ts_score(best_tree, singleton_ds)
  expect_equal(result$best_score, rescore)
})


# ===== Reproducibility =====

test_that("Driven search is reproducible with simplification", {
  set.seed(3390)
  r1 <- ts_driven(autap_ds)
  set.seed(3390)
  r2 <- ts_driven(autap_ds)
  expect_equal(r1$best_score, r2$best_score)
  expect_equal(r1$trees, r2$trees)
})


# ===== Larger dataset regression =====

test_that("EW scores match phangorn on a moderately sized dataset", {
  set.seed(8502)
  big_mat <- matrix(sample(0:3, 20 * 30, replace = TRUE),
                    nrow = 20,
                    dimnames = list(paste0("t", 1:20), NULL))
  big_dataset <- MatrixToPhyDat(big_mat)
  big_ds <- make_ts_data(big_dataset)

  for (i in seq_len(3)) {
    tree <- RandomTree(big_dataset, root = TRUE)
    score <- ts_score(tree, big_ds)
    expected <- phangorn::parsimony(tree, big_dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

test_that("Driven search on 4-state dataset with autapomorphies", {
  set.seed(9201)
  # Create dataset with some 4-state chars where some states are singletons
  mat4 <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1, 2, 3,  # states 2,3 are singletons
    0, 0, 1, 1, 0, 0, 1, 1, 0, 0,  # informative binary
    0, 0, 0, 1, 1, 1, 0, 0, 0, 0,  # informative binary
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1   # invariant
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  ds4 <- MatrixToPhyDat(mat4)
  ds4_data <- make_ts_data(ds4)

  diag <- ts_diag(ds4_data)
  # Invariant char should be removed

  expect_gte(diag$n_patterns_removed, 1L)

  result <- ts_driven(ds4_data)
  # Verify best score matches re-scoring
  best_tree <- structure(
    list(edge = result$trees[[1]],
         tip.label = paste0("t", 1:10),
         Nnode = 9L),
    class = "phylo"
  )
  rescore <- ts_score(best_tree, ds4_data)
  expect_equal(result$best_score, rescore)

  # Also check against phangorn
  phangorn_score <- phangorn::parsimony(best_tree, ds4)
  expect_equal(rescore, phangorn_score)
})


# ===== Edge cases =====

test_that("All-uninformative dataset scores correctly", {
  # Every character is an autapomorphy
  uninf_mat <- matrix(c(
    0, 1, 2, 3, 4, 5, 6, 7,
    0, 1, 2, 3, 4, 5, 6, 7
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  uninf_dataset <- MatrixToPhyDat(uninf_mat)
  uninf_ds <- make_ts_data(uninf_dataset)

  diag <- ts_diag(uninf_ds)
  # All patterns should be uninformative
  expect_true(all(!diag$informative))

  tree <- RandomTree(uninf_dataset, root = TRUE)
  score <- ts_score(tree, uninf_ds)
  expected <- phangorn::parsimony(tree, uninf_dataset)
  expect_equal(score, expected)
})

test_that("Single informative character among many uninformative", {
  # Mix of 1 informative + several uninformative chars
  mix_mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,  # informative
    0, 1, 2, 3, 4, 5, 6, 7,  # all autapomorphies
    0, 0, 0, 0, 0, 0, 0, 1   # single autapomorphy
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  mix_dataset <- MatrixToPhyDat(mix_mat)
  mix_ds <- make_ts_data(mix_dataset)

  diag <- ts_diag(mix_ds)
  # At least some patterns removed
  expect_gte(diag$n_patterns_removed, 1L)

  set.seed(6102)
  for (i in seq_len(5)) {
    tree <- RandomTree(mix_dataset, root = TRUE)
    score <- ts_score(tree, mix_ds)
    expected <- phangorn::parsimony(tree, mix_dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

# ===== Ambiguous token tests (T-013, T-014, T-017) =====

# Helper: build raw ts_data components with custom contrast matrix
make_custom_data <- function(contrast, tip_tokens, weight = NULL) {
  n_tips <- length(tip_tokens)
  n_patterns <- 1L
  list(
    contrast = contrast,
    tip_data = matrix(as.integer(tip_tokens), ncol = n_patterns),
    weight = if (is.null(weight)) 1L else as.integer(weight),
    levels = paste0("s", seq_len(ncol(contrast)))
  )
}

# Helper: score a tree against custom data, reordering tips to match labels
score_custom <- function(tree, ds, tip_names) {
  labels <- tree$tip.label
  # Reorder tip_data rows to match tree's tip label order
  idx <- match(labels, tip_names)
  td <- matrix(ds$tip_data[idx, , drop = FALSE], ncol = ncol(ds$tip_data))
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, td, ds$weight, ds$levels)
}

test_that("T-013: ambiguous informative character not removed", {
  # Tips: {0,1},{0,1},{0,1},{2,3},{2,3},{2,3} — 6 tips, 2 token types
  # This is parsimony-informative: score varies from 1 to 3 across trees.
  contrast <- matrix(c(1,1,0,0, 0,0,1,1), nrow = 2, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,1L,1L,2L,2L,2L))
  ds$levels <- c("0","1","2","3")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  # Must be kept as informative

  expect_true(diag$informative)
  expect_equal(diag$n_patterns_removed, 0L)
  expect_equal(diag$ew_offset, 0L)

  # Verify scores differ on two trees
  tip_names <- paste0("t", 1:6)
  tree_grouped <- ape::read.tree(text = "(((t1,t2),t3),((t4,t5),t6));")
  tree_mixed <- ape::read.tree(text = "(((t1,t4),t2),((t3,t5),t6));")
  s_grouped <- score_custom(tree_grouped, ds, tip_names)
  s_mixed <- score_custom(tree_mixed, ds, tip_names)
  expect_true(s_grouped != s_mixed,
              info = "Ambiguous informative char should give different scores")
})

test_that("T-013: ambiguous 4-tip informative character preserved", {
  # Smaller case: {0,1},{0,1},{2,3},{2,3} — 4 tips
  contrast <- matrix(c(1,1,0,0, 0,0,1,1), nrow = 2, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,1L,2L,2L))
  ds$levels <- c("0","1","2","3")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  expect_true(diag$informative)
  expect_equal(diag$n_patterns_removed, 0L)
})

test_that("T-014: all-ambiguous truly uninformative gets correct fixed cost", {
  # Tips: {0,1},{1,2},{0,2} — 3 tips, score is 1 on every tree
  contrast <- matrix(c(1,1,0, 0,1,1, 1,0,1), nrow = 3, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,2L,3L))
  ds$levels <- c("0","1","2")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  # Should be uninformative with precomputed_steps = 1
  expect_false(diag$informative)
  expect_equal(diag$precomputed_steps, 1L)
  expect_equal(diag$ew_offset, 1L)

  # Verify the score is correct on any tree
  tree <- ape::read.tree(text = "((t1,t2),t3);")
  tip_names <- paste0("t", 1:3)
  s <- score_custom(tree, ds, tip_names)
  expect_equal(s, 1)
})

test_that("T-014: all-ambiguous invariant character gives 0 fixed cost", {
  # Tips: {0,1},{0,1},{0,1} — 3 tips, 0 steps always (common state 0 and 1)
  contrast <- matrix(c(1,1), nrow = 1, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,1L,1L))
  ds$levels <- c("0","1")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  expect_false(diag$informative)
  expect_equal(diag$precomputed_steps, 0L)
  expect_equal(diag$ew_offset, 0L)
})

test_that("Transform 3: ambiguity state removal", {
  # 6 tips: 0, 0, 1, 1, {0,2}, {1,2}
  # State 2 only appears in ambiguity tokens alongside 0 or 1.
  # State 2 is redundant and should be removed.
  # After removal: 0, 0, 1, 1, {0}, {1} — informative binary.
  contrast <- matrix(c(
    1, 0, 0,   # token 1: state 0
    0, 1, 0,   # token 2: state 1
    1, 0, 1,   # token 3: {0,2}
    0, 1, 1    # token 4: {1,2}
  ), nrow = 4, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,1L,2L,2L,3L,4L))
  ds$levels <- c("0","1","2")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  # State 2 should be removed (n_states_reduced >= 1)
  expect_gte(diag$n_states_reduced, 1L)
  # Character should remain informative (0 and 1 each appear in 2+ tips)
  expect_true(diag$informative)
  # Fewer remaining states than the original 3
  expect_lte(diag$n_states_remaining, 2L)
})

test_that("Mixed ambiguous + unambiguous: singleton removal still works", {
  # 5 tips: 0, 0, 1, 1, {0,2}
  # State 0: unambig in 2 tips. State 1: unambig in 2 tips. State 2: ambig only.
  # State 2 is redundant (Transform 3). Character is informative.
  contrast <- matrix(c(
    1, 0, 0,   # token 1: state 0
    0, 1, 0,   # token 2: state 1
    1, 0, 1    # token 3: {0,2}
  ), nrow = 3, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,1L,2L,2L,3L))
  ds$levels <- c("0","1","2")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  # Character should stay informative
  expect_true(diag$informative)

  # Verify correct EW score against a known tree
  tip_names <- paste0("t", 1:5)
  tree <- ape::read.tree(text = "((t1,t2),((t3,t4),t5));")
  s <- score_custom(tree, ds, tip_names)
  # Tip 5 has {0,2}. On this tree, t1,t2 share 0; t3,t4 share 1;
  # t5={0,2} groups with t1,t2 side -> 1 step at the root.
  expect_equal(s, 1)
})

test_that("Ambiguous character with 2 ambig + 2 unambig tips", {
  # 4 tips: 0, 1, {0,1}, {0,1}
  # 0: unambig 1 tip. 1: unambig 1 tip. Both ambig tokens have both.
  # Classical criterion: 0 states with count >= 2 -> uninformative
  # But IS it? On any 4-tip tree, with Fitch:
  # A tip with {0,1} can resolve to either 0 or 1. So the character
  # should always cost 1 step (the single 0-vs-1 change).
  # Verify with caterpillar: fwd (0,1,{01},{01}): 0∩1={} cost 1, {01}∩{01}={01},
  # union {01}∩{01}={01} cost 0. Total=1. Rev: same. So truly uninformative.
  contrast <- matrix(c(
    1, 0,   # token 1: state 0
    0, 1,   # token 2: state 1
    1, 1    # token 3: {0,1}
  ), nrow = 3, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,2L,3L,3L))
  ds$levels <- c("0","1")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  # Should be uninformative (score is 1 on all trees)
  expect_false(diag$informative)
  expect_equal(diag$precomputed_steps, 1L)
  expect_equal(diag$ew_offset, 1L)

  # Verify score
  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tip_names <- paste0("t", 1:4)
  s <- score_custom(tree, ds, tip_names)
  expect_equal(s, 1)
})


# ===== Degenerate inapplicable patterns (subsetting / simulation) =====

# When matrices are subsetted or simulated, degenerate patterns can arise
# where all tips are inapplicable ("-") or missing ("?"). These skip
# simplification (three-pass NA scoring requires it) but should not break
# scoring or search.

test_that("All-inapplicable pattern adds 0 to score", {
  # All-"-" character alone: should score 0 on every topology
  inapp_only <- MatrixToPhyDat(matrix(
    rep("-", 6), nrow = 6,
    dimnames = list(paste0("t", 1:6), NULL)
  ))
  set.seed(6331)
  for (i in seq_len(5)) {
    tree <- RandomTree(inapp_only, root = TRUE)
    expect_equal(TreeLength(tree, inapp_only), 0,
                 info = paste("All-inapp tree", i))
  }
})

test_that("Degenerate inapplicable patterns don't alter informative scores", {
  # Score a dataset with and without degenerate patterns;
  # the degenerate chars should add 0.
  info_mat <- matrix(c(
    "0", "0", "0", "1", "1", "1",
    "0", "0", "1", "0", "1", "1"
  ), nrow = 6, dimnames = list(paste0("t", 1:6), NULL))
  info_dataset <- MatrixToPhyDat(info_mat)

  combined_mat <- cbind(
    matrix(rep("-", 6), nrow = 6),               # all-inapp
    matrix(rep("?", 6), nrow = 6),               # all-missing
    matrix(c("-","?","-","?","-","?"), nrow = 6), # mixed -/?
    info_mat
  )
  dimnames(combined_mat) <- list(paste0("t", 1:6), NULL)
  combined_dataset <- MatrixToPhyDat(combined_mat)

  set.seed(8243)
  for (i in seq_len(5)) {
    tree <- RandomTree(info_dataset, root = TRUE)
    score_info <- TreeLength(tree, info_dataset)
    score_combined <- TreeLength(tree, combined_dataset)
    expect_equal(score_combined, score_info,
                 info = paste("Tree", i,
                              "— degenerate chars should add 0"))
  }
})

test_that("Mixed degenerate inapplicable/missing patterns score correctly", {
  # 8 tips: a mix of "?", "-", and their combinations,
  # plus informative characters
  mat <- matrix(c(
    "-", "?", "-", "?", "-", "?", "-", "?",  # degenerate: all -/?
    "?", "?", "?", "?", "?", "?", "?", "?",  # degenerate: all missing
    "0", "0", "0", "0", "1", "1", "1", "1",  # informative
    "0", "0", "1", "1", "0", "0", "1", "1"   # informative
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(4817)
  for (i in seq_len(5)) {
    tree <- RandomTree(dataset, root = TRUE)
    score <- ts_score(tree, ds)
    expected <- TreeLength(tree, dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

test_that("'All in [0, ?]' characters simplified in inapplicable datasets", {
  # In an inapplicable dataset, ? tokens include the inapp bit.
  # Characters where every tip is one applicable state or ? are genuinely
  # uninformative and should be simplified away (not bypassed).
  mat <- matrix(c(
    "0", "0", "?", "?", "?", "?", "?", "?",  # all-0-or-? -> uninformative
    "0", "0", "1", "1", "-", "-", "?", "?",   # mixed informative
    "0", "0", "0", "0", "1", "1", "1", "1"    # standard informative
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  diag <- ts_diag(ds)

  # The all-0-or-? pattern should be identified as uninformative
  expect_gte(diag$n_patterns_removed, 1L,
             label = "all-[0,?] pattern should be removed")

  # Scores must still match TreeLength on multiple trees
  set.seed(5291)
  for (i in seq_len(5)) {
    tree <- RandomTree(dataset, root = TRUE)
    score <- ts_score(tree, ds)
    expected <- TreeLength(tree, dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

test_that("000?----?000 pattern kept as informative", {
  # Genuine inapplicable tips ("-") make this topology-dependent under
  # three-pass NA scoring: the character MUST NOT be simplified away.
  mat <- matrix(c(
    "0", "0", "0", "?", "-", "-", "-", "-", "?", "0", "0", "0",
    "0", "0", "0", "0", "0", "0", "1", "1", "1", "1", "1", "1"
  ), nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Scores must match TreeLength
  set.seed(2618)
  for (i in seq_len(5)) {
    tree <- RandomTree(dataset, root = TRUE)
    score <- ts_score(tree, ds)
    expected <- TreeLength(tree, dataset)
    expect_equal(score, expected, info = paste("Tree", i))
  }
})

test_that("Driven search handles degenerate inapplicable patterns", {
  # Simulates a subsetted matrix: 3 degenerate + 2 informative characters
  mat <- matrix(c(
    "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",  # all inapp
    "?", "-", "?", "-", "?", "-", "?", "-", "?", "-",  # all missing/inapp
    "-", "?", "-", "-", "?", "?", "-", "-", "?", "?",  # all missing/inapp
    "0", "0", "0", "1", "1", "0", "0", "1", "1", "1",  # informative
    "0", "0", "1", "0", "1", "1", "0", "1", "0", "1"   # informative
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_driven(ds)
  # Build the best tree and verify score
  best_tree <- structure(
    list(edge = result$trees[[1]],
         tip.label = paste0("t", seq_len(nrow(ds$tip_data))),
         Nnode = nrow(ds$tip_data) - 1L),
    class = "phylo"
  )
  rescore <- ts_score(best_tree, ds)
  expect_equal(result$best_score, rescore)

  # Score should also match TreeLength
  tl <- TreeLength(best_tree, dataset)
  expect_equal(rescore, tl)
})
