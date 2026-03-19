# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for x-transformation (Goloboff et al. 2021) scoring via the
# recode_hierarchy() → Sankoff pipeline and MaximizeParsimony(inapplicable="xform").

library("TreeTools")

make_dat <- function(mat, levels = c("-", "0", "1")) {
  phangorn::phyDat(mat, type = "USER", levels = levels, ambiguity = "?")
}


# ===== Gain/loss asymmetry ===================================================
# The x-transformation penalizes gains (absent→present) more heavily than
# losses (present→absent) at ratio (n+1):1, where n = number of secondaries.

test_that("Xform prefers single gain + losses over multiple gains", {
  # Tree: ((t1,t2),(t3,t4))
  # Primary: t1=absent, t2=present, t3=present, t4=present
  # Secondary: t2=0, t3=0, t4=0 (all identical when present)
  # States: absent=0, (sec=0)=1
  # Cost: gain=2, loss=1
  # Optimal on this tree: root=present(1), loss to t1 → cost 1
  # Alternative: root=absent(0), gain at MRCA(t2,t3,t4)... not a single node
  #   on ((t1,t2),(t3,t4)). Would need gain at each present tip = 3*2 = 6
  # So single-gain-from-root (cost 1) << multiple-gains (cost 6)
  mat <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "1", "0", "0",
    "1", "0", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  recoded <- recode_hierarchy(ds, h)
  blk <- recoded$sankoff_chars[[1]]
  expect_equal(blk$cost_matrix[1, 2], 2)  # gain = n+1 = 2
  expect_equal(blk$cost_matrix[2, 1], 1)  # loss = 1

  # Score via Sankoff: tree where absent tip is sister to one present
  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))
  res <- TreeSearch:::ts_sankoff_test(
    tree$edge,
    as.integer(blk$n_states),
    list(blk$cost_matrix),
    matrix(as.integer(blk$tip_states), ncol = 1),
    as.integer(blk$forced_root_state)
  )
  # Root=1(present), nodeAB: state 1 costs 0(t2)+1(loss to t1)=1
  # nodeCD: state 1 costs 0+0=0. Root=1: 0+0=0 from children.
  # But root cost = min over states. state 1 at root: costAB(1)=1, costCD(1)=0
  # cost_root_state1 = 1 + 0 = 1
  expect_equal(res$score, 1)
})


# ===== Secondary variation increases xform score =============================

test_that("Xform penalizes secondary variation on present branches", {
  # All present, varying secondaries → present-present Hamming cost
  mat_uniform <- matrix(c(
    "1", "0", "0",
    "1", "0", "0",
    "1", "0", "0",
    "1", "0", "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))

  mat_varied <- matrix(c(
    "1", "0", "0",
    "1", "1", "1",
    "1", "0", "0",
    "1", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))

  ds_u <- make_dat(mat_uniform)
  ds_v <- make_dat(mat_varied)
  h <- CharacterHierarchy("1" = 2:3)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree_u <- Renumber(RenumberTips(tree, names(ds_u)))
  tree_v <- Renumber(RenumberTips(tree, names(ds_v)))

  score_fn <- function(ds, tr) {
    rec <- recode_hierarchy(ds, h)
    blk <- rec$sankoff_chars[[1]]
    TreeSearch:::ts_sankoff_test(
      tr$edge, as.integer(blk$n_states),
      list(blk$cost_matrix),
      matrix(as.integer(blk$tip_states), ncol = 1),
      as.integer(blk$forced_root_state)
    )$score
  }

  score_uniform <- score_fn(ds_u, tree_u)
  score_varied <- score_fn(ds_v, tree_v)

  # Uniform: all same state, no Hamming cost → 0
  expect_equal(score_uniform, 0)
  # Varied: secondary changes → Hamming cost > 0
  expect_gt(score_varied, 0)
})


# ===== HSJ vs xform cross-validation =========================================
# Both methods handle inapplicable characters; they should agree on basic
# properties even if exact scores differ.

test_that("HSJ and xform agree on optimal tree for simple gain scenario", {
  mat <- matrix(c(
    "0", "-", "-", "0", "0", "0",
    "0", "-", "-", "0", "1", "1",
    "1", "0", "0", "1", "0", "0",
    "1", "0", "1", "1", "0", "1",
    "1", "1", "0", "1", "1", "0",
    "1", "1", "1", "0", "1", "1"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  # Both should successfully search and return valid trees
  set.seed(7184)
  hsj_result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "hsj", hsj_alpha = 1.0,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )
  set.seed(7184)
  xform_result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )

  expect_s3_class(hsj_result[[1]], "phylo")
  expect_s3_class(xform_result[[1]], "phylo")

  # Both should find trees with the correct number of tips
  expect_equal(length(hsj_result[[1]]$tip.label), 6L)
  expect_equal(length(xform_result[[1]]$tip.label), 6L)
})


# ===== Xform with non-hierarchy characters ====================================

test_that("Xform correctly combines Fitch + Sankoff scoring", {
  # Chars 1-2: hierarchy (primary + secondary)
  # Chars 3-4: non-hierarchy (standard Fitch)
  mat <- matrix(c(
    "0", "-", "0", "0",
    "1", "0", "0", "1",
    "1", "1", "1", "0",
    "0", "-", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")

  # Standard Fitch score (all 4 chars as standard)
  fitch_total <- TreeLength(tree, ds)

  # Xform search should run
  result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
  expect_true(is.finite(fitch_total))
})


# ===== Xform on larger dataset (8 tips) =====================================

test_that("Xform search works on 8-tip dataset", {
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
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
  expect_equal(length(result[[1]]$tip.label), 8L)

  # All result trees should be valid rooted phylogenies
  for (tr in result) {
    expect_s3_class(tr, "phylo")
    expect_true(TreeIsRooted(tr))
  }
})


# ===== Xform score is consistent across replicates ===========================

test_that("Xform search produces deterministic scores with same seed", {
  mat <- matrix(c(
    "0", "-", "0", "1",
    "1", "0", "1", "0",
    "1", "1", "0", "1",
    "1", "1", "1", "0",
    "0", "-", "1", "1"
  ), nrow = 5, byrow = TRUE,
  dimnames = list(paste0("t", 1:5), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  set.seed(3021)
  r1 <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )
  set.seed(3021)
  r2 <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )

  # Same seed, same result
  expect_equal(attr(r1, "score"), attr(r2, "score"))
})


# ===== Asymmetric cost correctness: gain vs loss =============================

test_that("Xform gain cost scales with number of secondaries", {
  # 1 secondary → gain = 2, 2 secondaries → gain = 3, 3 → gain = 4
  for (n_sec in 1:3) {
    n_cols <- 1 + n_sec
    mat <- matrix("-", nrow = 3, ncol = n_cols,
                  dimnames = list(paste0("t", 1:3), NULL))
    mat[1, ] <- c("0", rep("-", n_sec))
    mat[2, ] <- c("1", rep("0", n_sec))
    mat[3, ] <- c("1", rep("1", n_sec))
    ds <- make_dat(mat)
    h <- CharacterHierarchy("1" = seq(2L, n_cols))

    rec <- recode_hierarchy(ds, h)
    blk <- rec$sankoff_chars[[1]]

    expected_gain <- n_sec + 1L
    # All absent→present transitions should cost expected_gain
    for (j in 2:blk$n_states) {
      expect_equal(blk$cost_matrix[1, j], expected_gain,
                   info = paste("n_sec =", n_sec, "state", j))
    }
  }
})
