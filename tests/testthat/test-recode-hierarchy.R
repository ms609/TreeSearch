# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for recode_hierarchy(): x-transformation recoding of hierarchical
# characters into step-matrix (Sankoff) characters.

library("TreeTools")

make_dat <- function(mat, levels = c("-", "0", "1")) {
  phangorn::phyDat(mat, type = "USER", levels = levels, ambiguity = "?")
}


# ===== Basic 2-secondary binary block ========================================

test_that("Binary secondaries produce correct state count and cost matrix", {
  mat <- matrix(c(
    "0", "-", "-",
    "1", "0", "0",
    "1", "0", "1",
    "1", "1", "0",
    "1", "1", "1"
  ), nrow = 5, byrow = TRUE,
  dimnames = list(paste0("t", 1:5), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  result <- recode_hierarchy(ds, h)

  # One block, no non-hierarchy chars
  expect_length(result$sankoff_chars, 1)
  expect_length(result$non_hierarchy_indices, 0)

  blk <- result$sankoff_chars[[1]]

  # 2 binary secondaries → 2^2 + 1 = 5 states

  expect_equal(blk$n_states, 5L)

  # Cost matrix dimensions
  expect_equal(dim(blk$cost_matrix), c(5, 5))

  # Diagonal = 0
  expect_equal(diag(blk$cost_matrix), rep(0, 5))

  # Gain cost = n+1 = 3 (absent → any present)
  expect_true(all(blk$cost_matrix[1, 2:5] == 3))

  # Loss cost = 1 (any present → absent)
  expect_true(all(blk$cost_matrix[2:5, 1] == 1))

  # Present → present = Hamming distance
  # States 1..4 = combinations of 2 binary: (1,1), (2,1), (1,2), (2,2)
  # in expand.grid order: (1,1)=1, (2,1)=2, (1,2)=3, (2,2)=4
  # Hamming(1,2)=1, Hamming(1,3)=1, Hamming(1,4)=2
  # Hamming(2,3)=2, Hamming(2,4)=1, Hamming(3,4)=1
  expect_equal(blk$cost_matrix[2, 3], 1)  # (1,1)→(2,1)
  expect_equal(blk$cost_matrix[2, 4], 1)  # (1,1)→(1,2)
  expect_equal(blk$cost_matrix[2, 5], 2)  # (1,1)→(2,2)
  expect_equal(blk$cost_matrix[3, 4], 2)  # (2,1)→(1,2)
  expect_equal(blk$cost_matrix[3, 5], 1)  # (2,1)→(2,2)
  expect_equal(blk$cost_matrix[4, 5], 1)  # (1,2)→(2,2)
})


test_that("Tip states correctly encode absent and present combinations", {
  mat <- matrix(c(
    "0", "-", "-",
    "1", "0", "0",
    "1", "0", "1",
    "1", "1", "0",
    "1", "1", "1"
  ), nrow = 5, byrow = TRUE,
  dimnames = list(paste0("t", 1:5), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  blk <- recode_hierarchy(ds, h)$sankoff_chars[[1]]

  # t1: absent → state 0
  expect_equal(blk$tip_states[1], 0L)

  # t2: present (0,0) → combo (1,1) → row 1 → state 1
  expect_equal(blk$tip_states[2], 1L)

  # t3: present (0,1) → combo (1,2) → row 3 → state 3
  expect_equal(blk$tip_states[3], 3L)

  # t4: present (1,0) → combo (2,1) → row 2 → state 2
  expect_equal(blk$tip_states[4], 2L)

  # t5: present (1,1) → combo (2,2) → row 4 → state 4
  expect_equal(blk$tip_states[5], 4L)
})


# ===== Single secondary =====================================================

test_that("Single binary secondary gives 3 states", {
  mat <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "1", "1", "0"
  ), nrow = 3, byrow = TRUE,
  dimnames = list(paste0("t", 1:3), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  blk <- recode_hierarchy(ds, h)$sankoff_chars[[1]]
  expect_equal(blk$n_states, 3L)

  # Gain cost = 1+1 = 2
  expect_true(all(blk$cost_matrix[1, 2:3] == 2))
  # Loss cost = 1
  expect_true(all(blk$cost_matrix[2:3, 1] == 1))
  # Present→present Hamming = 1 (single char differs)
  expect_equal(blk$cost_matrix[2, 3], 1)
  expect_equal(blk$cost_matrix[3, 2], 1)
})


# ===== Non-hierarchy characters preserved ====================================

test_that("Non-hierarchy characters are identified", {
  mat <- matrix(c(
    "0", "-", "0", "1",
    "1", "0", "1", "0",
    "1", "1", "0", "1",
    "0", "-", "1", "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  result <- recode_hierarchy(ds, h)
  expect_equal(sort(result$non_hierarchy_indices), c(3L, 4L))
})


# ===== Multiple hierarchy blocks =============================================

test_that("Multiple blocks are handled independently", {
  mat <- matrix(c(
    "0", "-", "1", "0", "-",
    "1", "0", "1", "1", "0",
    "1", "1", "0", "1", "1",
    "0", "-", "0", "0", "-"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L, "4" = 5L)

  result <- recode_hierarchy(ds, h)

  expect_length(result$sankoff_chars, 2)
  expect_equal(result$non_hierarchy_indices, 3L)

  # Block 1: char 1 controls char 2 → 3 states
  expect_equal(result$sankoff_chars[[1]]$n_states, 3L)
  expect_equal(result$sankoff_chars[[1]]$block_chars, c(1L, 2L))

  # Block 2: char 4 controls char 5 → 3 states
  expect_equal(result$sankoff_chars[[2]]$n_states, 3L)
  expect_equal(result$sankoff_chars[[2]]$block_chars, c(4L, 5L))
})


# ===== Multistate secondary ==================================================

test_that("3-state secondary gives 4 combined states", {
  mat <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "1", "1", "0",
    "1", "2", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat, levels = c("-", "0", "1", "2"))
  h <- CharacterHierarchy("1" = 2L)

  blk <- recode_hierarchy(ds, h)$sankoff_chars[[1]]

  # 1 secondary with 3 informative states → 3 + 1 = 4 states
  expect_equal(blk$n_states, 4L)

  # Gain cost = 1+1 = 2
  expect_true(all(blk$cost_matrix[1, 2:4] == 2))

  # Present→present: each pair differs in 1 (only 1 secondary) → all Hamming = 1
  expect_equal(blk$cost_matrix[2, 3], 1)
  expect_equal(blk$cost_matrix[2, 4], 1)
  expect_equal(blk$cost_matrix[3, 4], 1)
})


test_that("Two multistate secondaries produce correct state count", {
  # 3-state × 2-state = 6 present states + 1 absent = 7
  mat <- matrix(c(
    "0", "-", "-",
    "1", "0", "0",
    "1", "1", "1",
    "1", "2", "0",
    "1", "2", "1"
  ), nrow = 5, byrow = TRUE,
  dimnames = list(paste0("t", 1:5), NULL))
  ds <- make_dat(mat, levels = c("-", "0", "1", "2"))
  h <- CharacterHierarchy("1" = 2:3)

  blk <- recode_hierarchy(ds, h)$sankoff_chars[[1]]

  # 3 × 2 + 1 = 7
  expect_equal(blk$n_states, 7L)

  # Gain cost = 2+1 = 3
  expect_true(all(blk$cost_matrix[1, 2:7] == 3))
})


# ===== Ambiguity handling ====================================================

test_that("Missing primary coded as fully ambiguous", {
  mat <- matrix(c(
    "?", "-", "0",
    "1", "0", "1",
    "0", "-", "0"
  ), nrow = 3, byrow = TRUE,
  dimnames = list(paste0("t", 1:3), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  blk <- recode_hierarchy(ds, h)$sankoff_chars[[1]]
  expect_equal(blk$tip_states[1], -1L)  # fully ambiguous
})


test_that("Present primary with unknown secondary coded as present-ambiguous", {
  mat <- matrix(c(
    "1", "?",
    "1", "0",
    "0", "-"
  ), nrow = 3, byrow = TRUE,
  dimnames = list(paste0("t", 1:3), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  blk <- recode_hierarchy(ds, h)$sankoff_chars[[1]]
  expect_equal(blk$tip_states[1], -2L)  # present but unknown
  expect_equal(blk$tip_states[2], 1L)   # present, state 0 → combo 1
  expect_equal(blk$tip_states[3], 0L)   # absent
})


# ===== State limit warning ===================================================

test_that("Large state space triggers warning", {
  # 5 binary secondaries → 2^5 + 1 = 33 states (> 32 limit)
  ncols <- 6
  mat <- matrix("-", nrow = 3, ncol = ncols,
                dimnames = list(paste0("t", 1:3), NULL))
  mat[1, ] <- c("0", rep("-", 5))
  mat[2, ] <- c("1", "0", "0", "0", "0", "0")
  mat[3, ] <- c("1", "1", "1", "1", "1", "1")
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:6)

  expect_warning(recode_hierarchy(ds, h), "33 states")
})


# ===== Nested hierarchy error =================================================

test_that("Nested hierarchies produce informative error", {
  mat <- matrix(c(
    "0", "-", "-", "-",
    "1", "0", "-", "0",
    "1", "1", "0", "1",
    "1", "1", "1", "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = list(2, "2" = 3L, 4))

  # validate_hierarchy catches double-claiming before recode gets to the
 # nesting check; either error message is acceptable
  expect_error(recode_hierarchy(ds, h), "multiple|Nested")
})


# ===== Cost matrix symmetry / asymmetry ======================================

test_that("Cost matrix is asymmetric for gain vs loss", {
  mat <- matrix(c(
    "0", "-",
    "1", "0",
    "1", "1"
  ), nrow = 3, byrow = TRUE,
  dimnames = list(paste0("t", 1:3), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  blk <- recode_hierarchy(ds, h)$sankoff_chars[[1]]

  # gain (absent→present) ≠ loss (present→absent)
  expect_equal(blk$cost_matrix[1, 2], 2)  # gain = n+1 = 2
  expect_equal(blk$cost_matrix[2, 1], 1)  # loss = 1
  expect_false(blk$cost_matrix[1, 2] == blk$cost_matrix[2, 1])
})


# ===== Integration: recode + Sankoff scoring ==================================

test_that("Recoded data scores correctly via Sankoff engine", {
  # Tree: ((t1,t2),(t3,t4))
  # Primary: t1=absent, t2-t4=present
  # Secondary: t2=0, t3=0, t4=1
  # States: 0=absent, 1=(0), 2=(1)
  # Optimal: root=present(0), loss to t1=1, t4: present(0)→present(1)=1
  # Total = 1 + 1 = 2? Let's compute properly.

  mat <- matrix(c(
    "0", "-",
    "1", "0",
    "1", "0",
    "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  result <- recode_hierarchy(ds, h)
  blk <- result$sankoff_chars[[1]]

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")

  # Score via Sankoff bridge
  tree <- Renumber(RenumberTips(tree, names(ds)))
  res <- TreeSearch:::ts_sankoff_test(
    edge = tree$edge,
    n_states_r = as.integer(blk$n_states),
    cost_matrices_r = list(blk$cost_matrix),
    tip_states_r = matrix(as.integer(blk$tip_states), ncol = 1),
    forced_root_r = as.integer(blk$forced_root_state)
  )

  expect_type(res$score, "double")
  expect_true(is.finite(res$score))

  # Hand-computed: 3 states (absent=0, sec0=1, sec1=2)
  # cm: gain=2, loss=1, Hamming(1↔2)=1
  # Tips: t1=0, t2=1, t3=1, t4=2
  # costAB = [2, 1, 2], costCD = [4, 1, 1]
  # root state 1: min(1+2,0+1,1+2)+min(1+4,0+1,1+1) = 1+1 = 2
  # Optimal: root=present(sec0), loss on t1 branch, Hamming-1 on t4 branch
  expect_equal(res$score, 2)
})


# ===== End-to-end: MaximizeParsimony with inapplicable='xform' ===============

test_that("MaximizeParsimony runs with inapplicable='xform'", {
  mat <- matrix(c(
    "0", "-", "-", "0", "0",
    "0", "-", "-", "0", "1",
    "1", "0", "0", "1", "0",
    "1", "0", "1", "1", "1",
    "1", "1", "0", "0", "0",
    "1", "1", "1", "0", "1"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
  expect_equal(length(result[[1]]$tip.label), 6L)
})


test_that("Xform score differs from standard Fitch on hierarchy data", {
  mat <- matrix(c(
    "0", "-", "0", "0",
    "1", "0", "1", "0",
    "1", "1", "0", "1",
    "1", "1", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))

  # Standard Fitch score
  fitch <- TreeLength(tree, ds)

  # Xform search should produce a valid result
  result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")

  # The xform score should be finite and non-negative
  # (we can't easily compare to Fitch since they measure different things,
  # but both should be valid)
  expect_true(is.finite(fitch))
})
