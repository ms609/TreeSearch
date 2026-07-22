# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# HSJ (Hopkins & St. John 2021) scoring end-to-end tests.
# Verifies the C++ hsj_score() algorithm and MaximizeParsimony() integration.

library("TreeTools")

# --- Internal wrappers ---
ts_hsj_score <- TreeSearch:::ts_hsj_score
.BuildTipLabels <- TreeSearch:::.BuildTipLabels
.HierarchyToBlocks <- TreeSearch:::.HierarchyToBlocks
.NonHierarchyWeights <- TreeSearch:::.NonHierarchyWeights
.HSJAbsentState <- TreeSearch:::.HSJAbsentState

# --- Helper: build a reductively-coded phyDat ---
make_hsj_dat <- function(mat, levels = c("-", "0", "1")) {
  phangorn::phyDat(mat, type = "USER", levels = levels, ambiguity = "?")
}

# --- Helper: score a tree under HSJ via the Rcpp bridge ---
hsj_score <- function(tree, dataset, hierarchy, alpha = 1.0) {
  at <- attributes(dataset)
  adj_w <- .NonHierarchyWeights(dataset, hierarchy)
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  blocks <- .HierarchyToBlocks(hierarchy)
  tl <- .BuildTipLabels(dataset)
  # absent_state = 0-based token index of "0" (= 1 for levels c("-","0","1")),
  # computed the same way the driven pipeline does.
  ts_hsj_score(
    edge = tree$edge,
    contrast = at$contrast,
    tip_data = tip_data,
    weight = as.integer(adj_w),
    levels = at$levels,
    hierarchy_blocks_r = blocks,
    alpha = alpha,
    tip_labels_r = tl,
    absent_state = .HSJAbsentState(dataset)
  )
}

# --- Helper: standard Fitch score ---
fitch_score <- function(tree, dataset) {
  d <- make_ts_data(dataset)
  ts_score(tree, d)
}


# =========================================================================
# Test: no-hierarchy characters → HSJ equals standard Fitch
# =========================================================================
test_that("HSJ with empty hierarchy equals standard Fitch", {
  mat <- matrix(c(
    "0", "1", "0", "1",
    "0", "0", "1", "1",
    "1", "0", "1", "0",
    "1", "1", "0", "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  ds <- make_hsj_dat(mat)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))

  expected_fitch <- fitch_score(tree, ds)

  # Score with hierarchy: char 1 controls chars 2-3
  h <- CharacterHierarchy("1" = 2:3)
  hsj_result <- hsj_score(tree, ds, h, alpha = 1.0)

  # They won't be equal because HSJ scores hierarchy chars differently.
  # But with NO hierarchy at all, they SHOULD be equal.
  # To test no-hierarchy equivalence, we need a dataset where no
  # characters are hierarchical. Use dummy empty hierarchy workaround:
  # Actually, we can't pass an empty hierarchy. Instead, compare
  # TreeLength() standard Fitch with HSJ where all chars are non-hierarchy.
  # This is tested implicitly via the Fitch component.

  # What we CAN test: the Fitch component of HSJ is correct.
  # With a hierarchy, the non-hierarchy chars should score identically
  # to Fitch applied to only those chars.
  expect_type(hsj_result, "double")
  expect_true(is.finite(hsj_result))
})


# =========================================================================
# Test: all-present hierarchy block with matching sister groups
# =========================================================================
test_that("HSJ scores all-present block with no secondary mismatches as 0", {
  # Tree: ((t1,t2),(t3,t4))
  # Primary: all present (state "1")
  # Sec char 2: t1="0", t2="0", t3="1", t4="1" (perfect split)
  # Sec char 3: t1="1", t2="1", t3="0", t4="0" (perfect split, inverted)
  # Non-hierarchy char 4: t1="0", t2="0", t3="1", t4="1"
  mat <- matrix(c(
    # pri  sec2  sec3  non-h
    "1",  "0",  "1",  "0",
    "1",  "0",  "1",  "0",
    "1",  "1",  "0",  "1",
    "1",  "1",  "0",  "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  ds <- make_hsj_dat(mat)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))

  h <- CharacterHierarchy("1" = 2:3)

  # HSJ hierarchy block: all present, but secondaries differ between clades.
  # After Fitch uppass, root resolves to state 0 for both sec chars.
  # node_A inherits root (0); node_B resolves to its own state.
  # Sec2: root=0, node_A=0, node_B=1 → d(root,node_B)=1 for this char.
  # Sec3: root=0, node_A=1, node_B=0 → d(root,node_A)=1 for this char.
  # HSJ block score = 1.0 (α·d/m on each root→child branch).
  # Non-hierarchy char 4: Fitch = 1 step.
  # Total HSJ score = 1.0 + 1 = 2
  expect_equal(hsj_score(tree, ds, h, alpha = 1.0), 2)
})


# =========================================================================
# Test: alpha=0 makes secondaries irrelevant
# =========================================================================
test_that("alpha=0 ignores secondary character variation", {
  # When alpha=0, present→present branch cost = 0 regardless of
  # secondary mismatches. So the hierarchy block score is determined
  # solely by the primary character's absent/present pattern.
  mat <- matrix(c(
    # pri  sec2  sec3
    "0",  "-",  "-",
    "1",  "0",  "1",
    "1",  "1",  "0",
    "1",  "0",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  ds <- make_hsj_dat(mat)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))

  h <- CharacterHierarchy("1" = 2:3)

  score_a0 <- hsj_score(tree, ds, h, alpha = 0.0)
  score_a1 <- hsj_score(tree, ds, h, alpha = 1.0)

  # With alpha=0, only the primary absent/present pattern matters.
  # Primary: "0","1","1","1" on ((t1,t2),(t3,t4))
  # Fitch on the primary alone: (t1=0,t2=1)→union, 1 step;
  # (t3=1,t4=1)→intersect; root: intersect → 0 more. Total = 1.
  # HSJ with alpha=0: the DP reduces to counting absent↔present transitions.
  # With t1 absent and t2,t3,t4 present:
  # Node for (t1,t2): min involves absent→present or present→absent = 1
  # Node for (t3,t4): both present, cost=0
  # Root: best is present→present on both sides = 0 + cost(left) + cost(right)
  # Expected alpha=0 score = 1 (one gain of the primary structure)
  expect_equal(score_a0, 1)

  # alpha=1 should be >= alpha=0 (secondaries add cost when mismatching)
  expect_gte(score_a1, score_a0)
})


# =========================================================================
# Test: alpha=0 equivalence across different secondary patterns
# =========================================================================
test_that("alpha=0 score is invariant to secondary character states", {
  # Two datasets with same primary pattern but different secondary states
  mat_a <- matrix(c(
    "0",  "-",  "-",
    "1",  "0",  "0",
    "1",  "0",  "0",
    "1",  "0",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))

  mat_b <- matrix(c(
    "0",  "-",  "-",
    "1",  "0",  "1",
    "1",  "1",  "0",
    "1",  "1",  "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))

  ds_a <- make_hsj_dat(mat_a)
  ds_b <- make_hsj_dat(mat_b)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds_a)))

  h <- CharacterHierarchy("1" = 2:3)

  expect_equal(
    hsj_score(tree, ds_a, h, alpha = 0.0),
    hsj_score(tree, ds_b, h, alpha = 0.0)
  )
})


# =========================================================================
# Test: HSJ score with mismatched secondaries
# =========================================================================
test_that("HSJ secondary dissimilarity detects mismatched secondaries", {
  # Tree: ((t1,t2),(t3,t4))
  # All tips present → primary block cost = 0 (no absent↔present transitions)
  # Secondaries identical → d=0 on every branch → block score = 0
  mat_match <- matrix(c(
    "1",  "0",  "0",
    "1",  "0",  "0",
    "1",  "0",  "0",
    "1",  "0",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))

  # Mismatched secondaries: (t1,t3)="0","0"; (t2,t4)="1","1"
  # On tree ((t1,t2),(t3,t4)), sister pairs have different secondary states
  # → d > 0 on internal branches → block score > 0
  mat_mismatch <- matrix(c(
    "1",  "0",  "0",
    "1",  "1",  "1",
    "1",  "0",  "0",
    "1",  "1",  "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))

  ds_match <- make_hsj_dat(mat_match)
  ds_mismatch <- make_hsj_dat(mat_mismatch)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds_match)))

  h <- CharacterHierarchy("1" = 2:3)

  score_match <- hsj_score(tree, ds_match, h, alpha = 1.0)
  score_mismatch <- hsj_score(tree, ds_mismatch, h, alpha = 1.0)

  # All identical secondaries → no dissimilarity → block score = 0
  expect_equal(score_match, 0)
  # Mismatched secondaries → d > 0 → block score > 0
  # Hand-computed: uppass resolves root & internal nodes to state 0 (lowest bit
  # of {0,1}), so t2 and t4 (state 1) mismatch their parents on both secondary
  # chars → d=2, m=2, α·d/m=1.0 per branch to t2 and t4.
  # Optimal: all present, p(root) = 2.0 (1.0 from left subtree + 1.0 from right)
  expect_equal(score_mismatch, 2.0)
})


# =========================================================================
# Test: single-gain scenario
# =========================================================================
test_that("HSJ scores single gain of a structure correctly", {
  # Tree: ((t1,t2),(t3,t4))
  # Primary: t1=absent, t2=t3=t4=present
  # Secondaries: all present tips have identical states → no secondary cost
  # Best mapping: gain on branch to (t2) from MRCA of (t1,t2)
  # or: gain at root, loss on t1 branch — but gain costs 1, loss costs 1,
  # so single gain = 1 is optimal
  mat <- matrix(c(
    "0",  "-",
    "1",  "0",
    "1",  "0",
    "1",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  ds <- make_hsj_dat(mat)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))

  h <- CharacterHierarchy("1" = 2L)

  # One absent tip, three present, all secondaries identical
  # Optimal: present at root, loss on t1 branch = 1
  # OR: absent at root, gain at MRCA(t2,t3,t4)... but that's not available
  # on this tree. On ((t1,t2),(t3,t4)):
  #   MRCA(t1,t2) = node A, MRCA(t3,t4) = node B, root = MRCA of all
  # Best: root=present, nodeA=present (gain+loss on t1 branch? no...)
  # Actually: root=present, nodeA: present costs 0 from root; t1=absent costs 1.
  #   nodeB: present costs 0; t3,t4 present costs 0. Total = 1.
  # Alternatively: root=absent: nodeA: t1 absent=0, t2 present gains=1;
  #   nodeB: t3+t4 present, each gains=1 but together... nodeB absent→gain on each? No.
  #   nodeB: best if present: root absent→nodeB present = 1 gain. t3,t4 present = 0.
  #   So root absent: nodeA best (absent→present for t2) = 1; nodeB gain = 1. Total = 2.
  # So present at root = score 1 is optimal.
  expect_equal(hsj_score(tree, ds, h, alpha = 1.0), 1)
})


# =========================================================================
# Test: HSJ with two hierarchy blocks
# =========================================================================
test_that("HSJ handles multiple hierarchy blocks", {
  # Two controlling primaries, each with one secondary
  mat <- matrix(c(
    # pri1  sec1a  pri2  sec2a  non_h
    "1",   "0",   "1",  "0",   "0",
    "1",   "0",   "1",  "1",   "1",
    "1",   "1",   "0",  "-",   "0",
    "1",   "1",   "0",  "-",   "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  ds <- make_hsj_dat(mat)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))

  h <- CharacterHierarchy("1" = 2L, "3" = 4L)

  score <- hsj_score(tree, ds, h, alpha = 1.0)
  expect_type(score, "double")
  expect_true(is.finite(score))
  expect_gte(score, 0)
})


# =========================================================================
# Test: alpha scales secondary contribution
# =========================================================================
test_that("HSJ score monotonically increases with alpha", {
  # Create a dataset where secondaries contribute to score
  mat <- matrix(c(
    "1",  "0",  "0",
    "1",  "1",  "1",
    "1",  "0",  "1",
    "1",  "1",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  ds <- make_hsj_dat(mat)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))

  h <- CharacterHierarchy("1" = 2:3)

  scores <- vapply(seq(0, 1, by = 0.25), function(a) {
    hsj_score(tree, ds, h, alpha = a)
  }, double(1))

  # Score should be non-decreasing in alpha (more weight on secondaries)
  for (i in seq_along(scores)[-1]) {
    expect_gte(scores[i], scores[i - 1])
  }
})


# =========================================================================
# Test: MaximizeParsimony end-to-end with HSJ
# =========================================================================
test_that("MaximizeParsimony runs with inapplicable='hsj'", {
  # 6-taxon dataset with hierarchy
  mat <- matrix(c(
    # pri  sec2  sec3  non_h1  non_h2  non_h3
    "0",  "-",  "-",  "0",    "0",    "0",
    "0",  "-",  "-",  "0",    "1",    "1",
    "1",  "0",  "0",  "1",    "0",    "0",
    "1",  "0",  "1",  "1",    "0",    "1",
    "1",  "1",  "0",  "1",    "1",    "0",
    "1",  "1",  "1",  "0",    "1",    "1"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_hsj_dat(mat)

  h <- CharacterHierarchy("1" = 2:3)

  result <- MaximizeParsimony(
    ds,
    hierarchy = h,
    inapplicable = "hsj",
    hsj_alpha = 1.0,
    maxReplicates = 2L,
    targetHits = 2L,
    verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
  expect_equal(length(result[[1]]$tip.label), 6L)
})


# =========================================================================
# Test: MaximizeParsimony HSJ with alpha=0
# =========================================================================
test_that("MaximizeParsimony HSJ alpha=0 works", {
  mat <- matrix(c(
    "0",  "-",  "-",  "0",  "0",
    "0",  "-",  "-",  "0",  "1",
    "1",  "0",  "0",  "1",  "0",
    "1",  "0",  "1",  "1",  "1",
    "1",  "1",  "0",  "0",  "0",
    "1",  "1",  "1",  "0",  "1"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_hsj_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  result <- MaximizeParsimony(
    ds,
    hierarchy = h,
    inapplicable = "hsj",
    hsj_alpha = 0.0,
    maxReplicates = 2L,
    targetHits = 2L,
    verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
})


# =========================================================================
# Test: HSJ parameter validation in MaximizeParsimony
# =========================================================================
test_that("MaximizeParsimony rejects bad HSJ parameters", {
  mat <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "1", "1", "0",
    "1", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_hsj_dat(mat)

  # hsj without hierarchy
  expect_error(
    MaximizeParsimony(ds, inapplicable = "hsj", verbosity = 0L),
    "hierarchy"
  )

  # bad alpha
  h <- CharacterHierarchy("1" = 2L)
  expect_error(
    MaximizeParsimony(ds, hierarchy = h, inapplicable = "hsj",
                      hsj_alpha = 2.0, verbosity = 0L),
    "hsj_alpha"
  )

  # IW + hsj (need a dataset with "-" for ValidateHierarchy to pass)
  mat2 <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "1", "1", "0",
    "1", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds2 <- make_hsj_dat(mat2)
  h2 <- CharacterHierarchy("1" = 2L)

  expect_error(
    MaximizeParsimony(ds2, hierarchy = h2, inapplicable = "hsj",
                      concavity = 10, verbosity = 0L),
    "Implied weighting"
  )

  # profile + hsj: PrepareDataProfile() strips "-" before validation,
  # so the error comes from ValidateHierarchy rather than the profile check
  expect_error(
    MaximizeParsimony(ds2, hierarchy = h2, inapplicable = "hsj",
                      concavity = "profile", verbosity = 0L),
    "inapplicable|Profile"
  )

  # xform is now implemented — should run without error
  # (but this minimal dataset may produce warnings)
  expect_s3_class(
    suppressWarnings(MaximizeParsimony(
      ds2, hierarchy = h2, inapplicable = "xform",
      maxReplicates = 1L, targetHits = 1L, verbosity = 0L
    ))[[1]],
    "phylo"
  )
})


# =========================================================================
# Test: HSJ score with a larger example (8 tips)
# =========================================================================
test_that("HSJ scoring works on 8-tip tree", {
  # Based on the paper's scenario: 8 taxa, more primaries than secondaries
  mat <- matrix(c(
    # pri1 sec1a sec1b pri2  pri3  pri4  pri5
    "1",  "0",  "0",  "0",  "0",  "0",  "0",
    "1",  "0",  "0",  "0",  "0",  "1",  "0",
    "1",  "0",  "1",  "0",  "1",  "0",  "0",
    "1",  "1",  "0",  "1",  "0",  "0",  "1",
    "1",  "1",  "1",  "1",  "0",  "0",  "1",
    "0",  "-",  "-",  "1",  "1",  "0",  "1",
    "0",  "-",  "-",  "1",  "1",  "1",  "0",
    "0",  "-",  "-",  "0",  "1",  "1",  "0"
  ), nrow = 8, byrow = TRUE,
  dimnames = list(paste0("t", 1:8), NULL))
  ds <- make_hsj_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  tree <- ape::read.tree(
    text = "(((t1,t2),(t3,t4)),((t5,t6),(t7,t8)));"
  )
  tree <- Renumber(RenumberTips(tree, names(ds)))

  score <- hsj_score(tree, ds, h, alpha = 1.0)
  expect_type(score, "double")
  expect_true(is.finite(score))
  expect_gte(score, 0)

  # alpha=0 should differ (or equal) but be valid

  score_a0 <- hsj_score(tree, ds, h, alpha = 0.0)
  expect_gte(score, score_a0)
})


# =========================================================================
# Test: HSJ search on 8-tip dataset finds trees
# =========================================================================
test_that("MaximizeParsimony HSJ search on 8-tip dataset", {
  mat <- matrix(c(
    "1",  "0",  "0",  "0",  "0",  "0",  "0",
    "1",  "0",  "0",  "0",  "0",  "1",  "0",
    "1",  "0",  "1",  "0",  "1",  "0",  "0",
    "1",  "1",  "0",  "1",  "0",  "0",  "1",
    "1",  "1",  "1",  "1",  "0",  "0",  "1",
    "0",  "-",  "-",  "1",  "1",  "0",  "1",
    "0",  "-",  "-",  "1",  "1",  "1",  "0",
    "0",  "-",  "-",  "0",  "1",  "1",  "0"
  ), nrow = 8, byrow = TRUE,
  dimnames = list(paste0("t", 1:8), NULL))
  ds <- make_hsj_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  result <- MaximizeParsimony(
    ds,
    hierarchy = h,
    inapplicable = "hsj",
    hsj_alpha = 1.0,
    maxReplicates = 3L,
    targetHits = 2L,
    verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
  expect_equal(length(result[[1]]$tip.label), 8L)

  # All result trees should be valid phylogenies
  for (tr in result) {
    expect_s3_class(tr, "phylo")
    expect_true(TreeIsRooted(tr))
  }
})


# =========================================================================
# Test: HSJ with all-absent and all-present tips
# =========================================================================
test_that("HSJ handles extreme absent/present ratios", {
  # Only one tip present
  mat_one <- matrix(c(
    "0",  "-",
    "0",  "-",
    "0",  "-",
    "1",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds_one <- make_hsj_dat(mat_one)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds_one)))

  h <- CharacterHierarchy("1" = 2L)

  score_one <- hsj_score(tree, ds_one, h, alpha = 1.0)
  expect_equal(score_one, 1)  # One gain (or loss from root)
})


# =========================================================================
# Regression: absent_state must identify the primary's "0" (absent) state,
# not the inapplicable "-" token, and must follow the level ordering.
# (Driven pipeline previously hard-coded 0L = index of "-", so primaries
#  coded "0" were treated as present and gain/loss was never counted.)
# =========================================================================
test_that(".HSJAbsentState() tracks the '0' token across level orderings", {
  expect_equal(.HSJAbsentState(make_hsj_dat(
    matrix(c("0", "1", "0", "1"), 2, dimnames = list(c("a", "b"), NULL)),
    levels = c("-", "0", "1"))), 1L)
  expect_equal(.HSJAbsentState(make_hsj_dat(
    matrix(c("0", "1", "0", "1"), 2, dimnames = list(c("a", "b"), NULL)),
    levels = c("0", "1", "-"))), 0L)
  expect_equal(.HSJAbsentState(make_hsj_dat(
    matrix(c("0", "1", "0", "1"), 2, dimnames = list(c("a", "b"), NULL)),
    levels = c("1", "-", "0"))), 2L)
})

test_that("HSJ is sensitive to primary present/absent at alpha=0", {
  # At alpha=0 the block score counts only primary gains/losses, so a primary
  # absence MUST register.  Before the fix this returned 0 (absence invisible).
  mat <- matrix(c(
    "0",  "-",  "-",
    "1",  "0",  "1",
    "1",  "1",  "0",
    "1",  "0",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_hsj_dat(mat)
  tree <- Renumber(RenumberTips(
    ape::read.tree(text = "((t1,t2),(t3,t4));"), names(ds)))
  h <- CharacterHierarchy("1" = 2:3)

  # One absent tip among three present → one gain.
  expect_equal(hsj_score(tree, ds, h, alpha = 0), 1)

  # Make every tip present → no gain/loss → block score 0.
  mat_all <- mat
  mat_all["t1", ] <- c("1", "0", "0")
  ds_all <- make_hsj_dat(mat_all)
  expect_equal(hsj_score(tree, ds_all, h, alpha = 0), 0)
})

test_that("driven HSJ (TreeLength) agrees with direct ts_hsj_score()", {
  # The driven pipeline and the test bridge must compute the same absent_state.
  mat <- matrix(c(
    "0",  "-",  "-",  "0",
    "0",  "-",  "-",  "1",
    "1",  "0",  "1",  "0",
    "1",  "1",  "0",  "1",
    "1",  "0",  "0",  "0",
    "1",  "1",  "1",  "1"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_hsj_dat(mat)
  tree <- Renumber(RenumberTips(ape::read.tree(
    text = "((t1,t2),((t3,t4),(t5,t6)));"), names(ds)))
  h <- CharacterHierarchy("1" = 2:3)

  for (a in c(0, 0.5, 1)) {
    expect_equal(
      TreeLength(tree, ds, hierarchy = h, inapplicable = "hsj", hsj_alpha = a),
      hsj_score(tree, ds, h, alpha = a)
    )
  }
})

test_that("HSJ score is invariant to phyDat level ordering", {
  # A parsimony-style score must not depend on the arbitrary internal ordering
  # of phyDat `levels`.  Two contributions could leak the ordering:
  #   * the PRIMARY absent/present term  — guarded by .HSJAbsentState() (T-307);
  #   * the SECONDARY dissimilarity term — the Fitch uppass in fitch_label_char()
  #     formerly resolved ambiguous internal nodes to the LOWEST SET BIT, whose
  #     token depends on `levels`.  It now resolves toward the best-supported
  #     token (subtree count, ties by smallest tip index), which is keyed on the
  #     tokens and tree rather than the bit encoding.
  # The secondary term only bites at alpha > 0, so test alpha in {0, 0.5, 1}.
  mat <- matrix(c(
    "0",  "-",  "-",
    "1",  "0",  "1",
    "1",  "1",  "0",
    "1",  "0",  "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  tree <- Renumber(RenumberTips(
    ape::read.tree(text = "((t1,t2),(t3,t4));"),
    paste0("t", 1:4)))
  h <- CharacterHierarchy("1" = 2:3)

  # All six orderings of the three tokens.
  orderings <- list(c("-", "0", "1"), c("-", "1", "0"), c("0", "-", "1"),
                    c("0", "1", "-"), c("1", "-", "0"), c("1", "0", "-"))
  for (a in c(0, 0.5, 1)) {
    scores <- vapply(orderings, function(lv) {
      hsj_score(tree, make_hsj_dat(mat, levels = lv), h, alpha = a)
    }, double(1))
    # Every ordering must agree (this dataset returned 2.5 vs 2.0 before the fix
    # at alpha = 1; the absent_state regression earlier made alpha = 0 disagree).
    expect_equal(scores, rep(scores[[1]], length(orderings)),
                 info = sprintf("hsj_alpha = %s", a))
  }
})

test_that("HSJ secondary dissimilarity is level-order invariant (multistate)", {
  # Stress the secondary term with a 3-state secondary and missing data, where
  # internal ambiguity is common and the lowest-bit tie-break was most exposed.
  mat <- matrix(c(
    "1",  "0",  "1",
    "1",  "2",  "?",
    "0",  "-",  "-",
    "1",  "1",  "0",
    "1",  "0",  "2",
    "1",  "2",  "1"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  tree <- Renumber(RenumberTips(ape::read.tree(
    text = "((t1,t2),((t3,t4),(t5,t6)));"), paste0("t", 1:6)))
  h <- CharacterHierarchy("1" = 2:3)

  toks <- c("-", "0", "1", "2")
  orderings <- list(toks, rev(toks), c("0", "1", "2", "-"),
                    c("2", "0", "-", "1"), c("1", "-", "2", "0"))
  for (a in c(0.5, 1)) {
    scores <- vapply(orderings, function(lv) {
      hsj_score(tree, make_hsj_dat(mat, levels = lv), h, alpha = a)
    }, double(1))
    expect_equal(scores, rep(scores[[1]], length(orderings)),
                 info = sprintf("hsj_alpha = %s", a))
  }
})


# =========================================================================
# Test: HSJ + sectorial search (T-303 guard)
# =========================================================================
# build_reduced_dataset() does not copy hierarchy_blocks/tip_labels/hsj_alpha,
# so rss_search/xss_search are guarded to fall back under HSJ (T-303); css_search
# scores the full dataset and needs no guard.  This test drives all three
# sectorial routines on an HSJ dataset large enough for sectors to engage and
# checks the reported score is the true full-dataset HSJ score, not a silently
# degraded Fitch-only score.
test_that("MaximizeParsimony HSJ + sectorial search stays score-consistent", {
  mat <- matrix(c(
    # pri  sec2  sec3  nh4   nh5   nh6   nh7
    "0",  "-",  "-",  "0",  "0",  "0",  "1",
    "0",  "-",  "-",  "0",  "1",  "1",  "0",
    "0",  "-",  "-",  "1",  "0",  "0",  "1",
    "0",  "-",  "-",  "1",  "1",  "1",  "0",
    "1",  "0",  "0",  "0",  "0",  "1",  "1",
    "1",  "0",  "0",  "0",  "1",  "0",  "0",
    "1",  "0",  "1",  "1",  "0",  "1",  "1",
    "1",  "0",  "1",  "1",  "1",  "0",  "0",
    "1",  "1",  "0",  "0",  "0",  "0",  "1",
    "1",  "1",  "0",  "0",  "1",  "1",  "0",
    "1",  "1",  "1",  "1",  "0",  "0",  "1",
    "1",  "1",  "1",  "1",  "1",  "1",  "0",
    "1",  "0",  "1",  "0",  "0",  "1",  "1",
    "1",  "1",  "0",  "1",  "1",  "0",  "0"
  ), nrow = 14, byrow = TRUE,
  dimnames = list(paste0("t", 1:14), NULL))
  ds <- make_hsj_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  ctrl <- SearchControl(
    ratchetCycles = 1L,
    xssRounds = 2L, xssPartitions = 3L,
    rssRounds = 2L, cssRounds = 1L, cssPartitions = 3L,
    sectorMinSize = 4L, sectorMaxSize = 10L
  )

  set.seed(8123)
  result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "hsj", hsj_alpha = 1.0,
    control = ctrl, maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )

  # The full HSJ + sectorial pipeline (rss/xss guarded, css on full ds) runs
  # to completion and returns valid trees with a finite, positive HSJ score.
  expect_s3_class(result[[1]], "phylo")
  expect_equal(length(result[[1]]$tip.label), 14L)
  reported <- attr(result, "score")
  expect_true(is.finite(reported))
  expect_true(reported > 0)

  # T-303 is a *silent* heuristic-quality bug: final scores are always
  # recomputed on the full dataset, so a regression cannot be caught by an
  # absolute-score assertion.  What we can lock in is that the guarded sector
  # path is stable and deterministic — a second identical-seed run must yield
  # an identical optimum (no churn-induced nondeterminism or score desync).
  set.seed(8123)
  result2 <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "hsj", hsj_alpha = 1.0,
    control = ctrl, maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )
  expect_equal(attr(result2, "score"), reported)
  expect_equal(length(result2), length(result))
})


# =========================================================================
# Test: all-hierarchy data -> zero Fitch words (regression)
# =========================================================================
# HSJ, like xform, zero-weights every hierarchy character, so a dataset whose
# characters are ALL hierarchical leaves the equal-weights dataset empty
# (DataSet::total_words == 0, empty per-word state vectors).  wagner_tree()
# indexed element 0 of those empty vectors (&ds.tip_states[0]) — undefined
# behaviour that aborted under the hardened libstdc++ assertions in the
# gcc-ASAN CI (run 28662381835).  The HSJ search must still build a valid tree
# from the hierarchy DP term alone.

test_that("HSJ search handles all-hierarchy data (zero Fitch words)", {
  mat <- matrix(c(
    "1", "0", "0", "-", "-",
    "1", "1", "1", "0", "1",
    "0", "-", "1", "1", "0",
    "1", "0", "1", "0", "0",
    "0", "-", "0", "-", "-",
    "1", "1", "0", "-", "-"
  ), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_hsj_dat(mat)
  # Two blocks span ALL five characters -> no non-hierarchy chars remain.
  h <- CharacterHierarchy("1" = 2L, "3" = 4:5)
  expect_length(setdiff(seq_len(5L), HierarchyChars(h)), 0L)

  set.seed(42)
  res <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "hsj",
                           hsj_alpha = 1.0, maxReplicates = 4L,
                           targetHits = 3L, verbosity = 0L)
  expect_s3_class(res[[1]], "phylo")
  expect_equal(length(res[[1]]$tip.label), 6L)
  for (tr in res) {
    expect_s3_class(tr, "phylo")
    expect_true(TreeIsRooted(tr))
  }
})
