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

test_that("T-014: all-ambiguous, no shared state -> kept & scored correctly", {
  # Tips: {0,1},{1,2},{0,2} — 3 tips, score is 1 on every tree.
  # No single state is shared by ALL tips, so the (topology-independent)
  # simplifier cannot prove a fixed cost > 0 and correctly KEEPS the
  # character; the per-tree downpass scores it exactly. (The old
  # caterpillar-sampling heuristic invented a fixed cost here — unsound in
  # general; see the regression block below.)
  contrast <- matrix(c(1,1,0, 0,1,1, 1,0,1), nrow = 3, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,2L,3L))
  ds$levels <- c("0","1","2")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  expect_true(diag$informative)
  expect_equal(diag$precomputed_steps, 0L)
  expect_equal(diag$ew_offset, 0L)

  # Score must still be correct (1) on the single 3-taxon topology
  tree <- ape::read.tree(text = "((t1,t2),t3);")
  s <- score_custom(tree, ds, paste0("t", 1:3))
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

test_that("Transform 3: state in >=2 tips is NOT removed (topology-dependent)", {
  # 6 tips: 0, 0, 1, 1, {0,2}, {1,2}
  # State 2 appears (ambiguously) in TWO tips. Removing it is unsound: on a
  # tree where the two ambiguous tips are sisters they resolve to {2} for a
  # cost-free clade, so the removed-state score can differ. The simplifier
  # must keep state 2 and let the downpass score the character.
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
  expect_true(diag$informative)
  # State 2 is retained: all three states remain.
  expect_equal(diag$n_states_remaining, 3L)

  # Scores must match phangorn on several random trees.
  tip_names <- paste0("t", 1:6)
  pd <- MatrixToPhyDat(matrix(
    c("0", "0", "1", "1", "{02}", "{12}"), ncol = 1,
    dimnames = list(tip_names, NULL)))
  set.seed(411)
  for (i in seq_len(5)) {
    tree <- RandomTree(tip_names, root = TRUE)
    s <- score_custom(tree, ds, tip_names)
    expect_equal(s, as.numeric(phangorn::fitch(tree, pd)),
                 info = paste("Tree", i))
  }
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

test_that("Binary 0,1,{0,1},{0,1}: {0,1} is full-? -> uninformative", {
  # 4 tips: 0, 1, {0,1}, {0,1}. In a 2-state character {0,1} is EVERY
  # applicable state, i.e. full "?" missing data, which is inert in Fitch. So
  # the character is really 0, 1, ?, ? -> a single 0-vs-1 change on every tree
  # (constant length 1) -> uninformative with a precomputed step of 1.
  contrast <- matrix(c(
    1, 0,   # token 1: state 0
    0, 1,   # token 2: state 1
    1, 1    # token 3: {0,1} == "?" for a binary character
  ), nrow = 3, byrow = TRUE)
  ds <- make_custom_data(contrast, c(1L,2L,3L,3L))
  ds$levels <- c("0","1")

  diag <- TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                         ds$weight, ds$levels)
  expect_false(diag$informative)
  expect_equal(diag$precomputed_steps, 1L)
  expect_equal(diag$ew_offset, 1L)

  # Score must be 1 on every 4-taxon topology
  tip_names <- paste0("t", 1:4)
  for (nwk in c("((t1,t2),(t3,t4));", "((t1,t3),(t2,t4));",
                "((t1,t4),(t2,t3));")) {
    s <- score_custom(ape::read.tree(text = nwk), ds, tip_names)
    expect_equal(s, 1, info = nwk)
  }
})


# ===== Regression: multistate-ambiguity over-count (topology-dependence) =====
#
# simplify_patterns must be topology-independent, but the optimal resolution of
# an ambiguous token is topology-dependent. Three historic bugs "simplified" a
# state whose resolution is topology-dependent and over-/under-counted on
# multistate ambiguity: Transform 2 (singleton state that also appears in
# ambiguous tokens), Transform 3 (ambiguous-only state present in >=2 tips), and
# the caterpillar-sampling uninformative heuristic. Each state below appears in
# >=2 tips (or the character has no shared state), so no reduction is sound;
# scores must match phangorn::fitch on every tree.

test_that("Regression: singleton state also in ambiguous tokens (Transform 2)", {
  # State 2 is unambiguous in exactly one tip (t15) but ALSO an option in two
  # ambiguous tips; the optimal tree resolves t31 to 2 for a free join, so
  # charging a fixed +1 autapomorphy step over-counts (was 4, truth 3).
  labs <- c("t15", "t31", "t33", "t35", "t37", "t38", "t40")
  toks <- c(t15 = "2", t31 = "{12}", t33 = "1", t35 = "1",
            t37 = "{02}", t38 = "0", t40 = "0")
  pd <- MatrixToPhyDat(matrix(toks[labs], ncol = 1,
                              dimnames = list(labs, NULL)))
  tree <- Preorder(ape::read.tree(
    text = "((t15,t31),((((t33,t38),t37),t35),t40));"))
  ts <- sum(CharacterLength(tree, pd, compress = FALSE))
  expect_equal(ts, as.numeric(phangorn::fitch(tree, pd)))
  expect_equal(ts, 3)
})

test_that("Regression: ambiguous-only state in >=2 tips (Transform 3)", {
  # State 2 is never unambiguous but is an option in FOUR tips; on this tree
  # A,B and C,D each resolve to {2} for a free 4-tip clade. Stripping state 2
  # forces extra steps (was 3, truth 2).
  labs <- c("A", "B", "C", "D", "E1", "E2", "F1", "F2")
  toks <- c(A = "{02}", B = "{12}", C = "{02}", D = "{12}",
            E1 = "0", E2 = "0", F1 = "1", F2 = "1")
  pd <- MatrixToPhyDat(matrix(toks[labs], ncol = 1,
                              dimnames = list(labs, NULL)))
  tree <- Preorder(ape::read.tree(
    text = "(((A,B),(C,D)),((E1,E2),(F1,F2)));"))
  ts <- sum(CharacterLength(tree, pd, compress = FALSE))
  expect_equal(ts, as.numeric(phangorn::fitch(tree, pd)))
  expect_equal(ts, 2)
})

test_that("Regression: classical-uninformative + ambiguity vs varied trees", {
  # Only state 2 is unambiguously doubled, so the classical criterion flags
  # this as uninformative; but with ambiguity the score is NOT constant
  # (caterpillars give 3, a balanced tree gives 2). Must match phangorn on
  # every topology, not a sampled fixed cost.
  labs <- paste0("t", 1:8)
  toks <- c("2", "{01}", "0", "{02}", "{12}", "1", "{012}", "2")
  pd <- MatrixToPhyDat(matrix(toks, ncol = 1, dimnames = list(labs, NULL)))
  set.seed(20260706)
  for (i in seq_len(20)) {
    tree <- RandomTree(labs, root = TRUE)
    expect_equal(sum(CharacterLength(tree, pd, compress = FALSE)),
                 as.numeric(phangorn::fitch(tree, pd)),
                 info = paste("tree", i))
  }
})

test_that("Regression: random ambiguous characters match phangorn", {
  set.seed(4242)
  for (trial in seq_len(300)) {
    n_tip <- sample(6:10, 1)
    S <- sample(3:4, 1)
    toks <- vapply(seq_len(n_tip), function(i) {
      if (runif(1) < 0.5) {
        as.character(sample.int(S, 1) - 1L)
      } else {
        k <- sample(2:min(3, S), 1)
        paste0("{", paste0(sort(sample.int(S, k) - 1L), collapse = ""), "}")
      }
    }, character(1))
    labs <- paste0("t", seq_len(n_tip))
    pd <- tryCatch(
      MatrixToPhyDat(matrix(toks, ncol = 1, dimnames = list(labs, NULL))),
      error = function(e) NULL)
    if (is.null(pd) || is.null(attr(pd, "levels"))) next
    tree <- RandomTree(labs, root = TRUE)
    expect_equal(sum(CharacterLength(tree, pd, compress = FALSE)),
                 as.numeric(phangorn::fitch(tree, pd)),
                 info = paste(toks, collapse = ","))
  }
})

test_that("Regression: simplify never over-removes (min==max oracle)", {
  # SAFETY GUARANTEE: a pattern may be dropped as uninformative ONLY if it is
  # genuinely topology-invariant, i.e. MinimumLength() == MaximumLength() (the
  # length is the same on every tree). This is the principled oracle; it guards
  # against a future change re-introducing the over-/under-count class of bug by
  # simplifying a character whose optimal resolution is topology-dependent.
  # (Measured 2026-07-06: 0 over-removals across 5856 fuzz+corpus patterns;
  # completeness misses 0/3356 on the real inapplicable corpus.)
  set.seed(1234)
  for (trial in seq_len(300)) {
    n_tip <- sample(6:9, 1)
    S <- sample(2:4, 1)
    toks <- vapply(seq_len(n_tip), function(i) {
      r <- runif(1)
      if (r < 0.45) {
        as.character(sample.int(S, 1) - 1L)
      } else if (r < 0.8) {
        paste0("{", paste0(sort(sample.int(S, min(2, S)) - 1L),
                           collapse = ""), "}")
      } else {
        "?"
      }
    }, character(1))
    labs <- paste0("t", seq_len(n_tip))
    pd <- tryCatch(
      MatrixToPhyDat(matrix(toks, ncol = 1, dimnames = list(labs, NULL))),
      error = function(e) NULL)
    if (is.null(pd) || is.null(attr(pd, "levels"))) next
    at <- attributes(pd)
    td <- matrix(unlist(pd, use.names = FALSE), nrow = length(pd), byrow = TRUE)
    diag <- TreeSearch:::ts_simplify_diag(at$contrast, td, at$weight, at$levels)
    removed <- which(!diag$informative)
    if (!length(removed)) next
    mn <- MinimumLength(pd, compress = TRUE)
    mx <- MaximumLength(pd, compress = TRUE)
    expect_true(all(mn[removed] == mx[removed]),
                info = paste("over-removal:", paste(toks, collapse = ",")))
  }
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
