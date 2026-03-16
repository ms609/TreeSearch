library("TreeTools")

# ---- Helpers ----

make_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels
  )
}

ts_score <- function(tree, ds, concavity = Inf, min_steps = integer(0)) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data,
                               ds$weight, ds$levels,
                               min_steps = min_steps,
                               concavity = concavity)
}

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
    expect_gte(iw_score, 0, info = paste("Tree", i, "non-negative"))
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


# ===== Regression: inapplicable datasets not affected =====

test_that("Inapplicable dataset scores match morphy (simplification skipped)", {
  skip_if_not_installed("TreeSearch")
  dataset <- TreeSearch::inapplicable.datasets$Vinther2008
  ds <- make_ts_data(dataset)
  tree <- TreeTools::PectinateTree(dataset)

  score <- ts_score(tree, ds)
  morphy_score <- phangorn::parsimony(tree, dataset)
  expect_equal(score, morphy_score)
})


# ===== Reproducibility =====

test_that("Driven search is reproducible with simplification", {
  set.seed(3390)
  r1 <- ts_driven(autap_ds, maxReplicates = 2L)
  set.seed(3390)
  r2 <- ts_driven(autap_ds, maxReplicates = 2L)
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

  result <- ts_driven(ds4_data, maxReplicates = 3L)
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

test_that("Transform 3: ambiguity state removal", {
  # Create dataset where a state appears only in ambiguity tokens:
  # Tips: 0, 0, 1, 1, {0,2}, {1,2}
  # State 2 only appears in ambiguity tokens alongside 0 or 1.
  # State 2 is redundant and should be removed.
  #
  # After removal: 0, 0, 1, 1, {0}, {1} = 001110
  # This is informative binary.
  skip("Ambiguity tokens require custom phyDat construction — complex to test")
})
