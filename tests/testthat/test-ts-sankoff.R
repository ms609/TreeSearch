# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Sankoff parsimony optimization engine unit tests.
# Hand-computed examples verify downpass scoring, root forcing, and uppass
# state reconstruction.

sankoff_test <- function(tree, n_states, cost_matrices, tip_states,
                         forced_root = rep(-1L, length(n_states))) {
  TreeSearch:::ts_sankoff_test(
    tree$edge,
    as.integer(n_states),
    cost_matrices,
    matrix(as.integer(tip_states), nrow = length(tree$tip.label)),
    as.integer(forced_root)
  )
}

# Helper: symmetric unit-cost matrix (Fitch-equivalent)
fitch_cost <- function(k) {
  m <- matrix(1, k, k)
  diag(m) <- 0
  m
}

# Helper: 4-tip balanced tree ((1,2),(3,4))
tree4 <- ape::read.tree(text = "((t1,t2),(t3,t4));")

# Helper: 5-tip pectinate tree (((((t1,t2),t3),t4),t5))
tree5 <- ape::read.tree(text = "(((t1,t2),t3),(t4,t5));")


# ===== Fitch equivalence (binary, symmetric unit cost) ====================

test_that("Sankoff matches Fitch for binary symmetric cost", {
  # A=0, B=0, C=1, D=1 on ((A,B),(C,D)) -> 1 step
  res <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(fitch_cost(2)),
    tip_states = c(0, 0, 1, 1)
  )
  expect_equal(res$score, 1)
  expect_equal(as.numeric(res$per_char), 1)
})

test_that("Sankoff matches Fitch for 3-state symmetric cost", {
  # A=0, B=1, C=2, D=0 on ((A,B),(C,D)) -> 2 steps
  res <- sankoff_test(
    tree4,
    n_states = 3L,
    cost_matrices = list(fitch_cost(3)),
    tip_states = c(0, 1, 2, 0)
  )
  expect_equal(res$score, 2)
})

test_that("No change needed -> score 0", {
  # All tips same state
  res <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(fitch_cost(2)),
    tip_states = c(0, 0, 0, 0)
  )
  expect_equal(res$score, 0)
})


# ===== Ordered (linear) cost matrix ======================================

test_that("Ordered character costs proportional to distance", {
  # 3-state ordered: cost[i][j] = |i - j|
  ordered3 <- matrix(c(0, 1, 2,
                        1, 0, 1,
                        2, 1, 0), 3, 3, byrow = TRUE)

  # A=0, B=0, C=2, D=2 -> optimal: AB=0, CD=2, root=1; cost = 0+0+1+1 = 2
  res <- sankoff_test(
    tree4,
    n_states = 3L,
    cost_matrices = list(ordered3),
    tip_states = c(0, 0, 2, 2)
  )
  expect_equal(res$score, 2)

  # A=0, B=2, C=0, D=2 -> optimal: internal nodes at state 1; cost = 4
  res2 <- sankoff_test(
    tree4,
    n_states = 3L,
    cost_matrices = list(ordered3),
    tip_states = c(0, 2, 0, 2)
  )
  expect_equal(res2$score, 4)
})


# ===== Asymmetric costs (gain/loss for x-transformation) ==================

test_that("Asymmetric gain:loss costs scored correctly", {
  # gain = 3, loss = 1 (from Goloboff x-transformation pattern)
  asym <- matrix(c(0, 3,
                    1, 0), 2, 2, byrow = TRUE)

  # A=0, B=0, C=1, D=1 -> root=1 (loss to AB): cost = 1; CD=1: cost = 0
  # Total = 1
  res <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(asym),
    tip_states = c(0, 0, 1, 1)
  )
  expect_equal(res$score, 1)

  # Reverse: A=1, B=1, C=0, D=0 -> root=1 (loss to CD): cost = 1
  res_rev <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(asym),
    tip_states = c(1, 1, 0, 0)
  )
  expect_equal(res_rev$score, 1)

  # A=0, B=1, C=0, D=0 -> optimal: root=1, AB=1, CD=0
  #   root→AB: 1→1=0, root→CD: 1→0=1(loss),
  #   AB→A: 1→0=1(loss), AB→B: 1→1=0, CD→C: 0→0=0, CD→D: 0→0=0
  #   Total = 0+1+1+0+0+0 = 2
  res3 <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(asym),
    tip_states = c(0, 1, 0, 0)
  )
  expect_equal(res3$score, 2)
})


# ===== Root forcing ======================================================

test_that("Forced root state changes score", {
  # gain = 3, loss = 1
  asym <- matrix(c(0, 3,
                    1, 0), 2, 2, byrow = TRUE)
  # A=0, B=0, C=1, D=1

  # Unconstrained: score = 1 (root=1, loss to AB)
  res_free <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(asym),
    tip_states = c(0, 0, 1, 1),
    forced_root = -1L
  )
  expect_equal(res_free$score, 1)

  # Force root = 0: root=0, AB=0 (cost 0), CD=1 (gain cost 3)
  res_r0 <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(asym),
    tip_states = c(0, 0, 1, 1),
    forced_root = 0L
  )
  expect_equal(res_r0$score, 3)

  # Force root = 1: root=1, AB=0 (loss cost 1), CD=1 (cost 0) -> 1
  res_r1 <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(asym),
    tip_states = c(0, 0, 1, 1),
    forced_root = 1L
  )
  expect_equal(res_r1$score, 1)
})


# ===== Uppass (state reconstruction) ======================================

test_that("Uppass assigns correct optimal states", {
  # Binary symmetric cost, A=0, B=0, C=1, D=1
  res <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(fitch_cost(2)),
    tip_states = c(0, 0, 1, 1)
  )

  # R node numbering: tips 1:4, root = 5, internal nodes 5,6,7
  # C++ numbering: tips 0:3, root = 4 (= n_tip), internal 4,5,6
  # optimal_states is n_node rows × n_chars cols (0-indexed states)
  opt <- res$optimal_states[, 1] # first (only) character

  # Tips must match observed

  expect_equal(opt[1], 0L)  # t1
  expect_equal(opt[2], 0L)  # t2
  expect_equal(opt[3], 1L)  # t3
  expect_equal(opt[4], 1L)  # t4

  # AB internal must be 0 (both children are 0)
  # CD internal must be 1 (both children are 1)
  # Root: either 0 or 1 (both cost the same); just check it's valid
  expect_true(opt[5] %in% c(0L, 1L))  # root
})

test_that("Uppass respects forced root state", {
  asym <- matrix(c(0, 3, 1, 0), 2, 2, byrow = TRUE)
  res <- sankoff_test(
    tree4,
    n_states = 2L,
    cost_matrices = list(asym),
    tip_states = c(0, 0, 1, 1),
    forced_root = 0L
  )

  opt <- res$optimal_states[, 1]
  expect_equal(opt[5], 0L)  # root forced to state 0
})


# ===== Multiple characters ================================================

test_that("Multi-character scoring sums correctly", {
  # Char 1: binary, A=0,B=0,C=1,D=1 -> 1 step
  # Char 2: binary, A=0,B=1,C=0,D=1 -> 2 steps
  res <- sankoff_test(
    tree4,
    n_states = c(2L, 2L),
    cost_matrices = list(fitch_cost(2), fitch_cost(2)),
    tip_states = cbind(c(0, 0, 1, 1), c(0, 1, 0, 1))
  )
  expect_equal(res$score, 3)
  expect_equal(as.numeric(res$per_char), c(1, 2))
})

test_that("Characters with different n_states handled correctly", {
  # Char 1: 2-state, cost 1; Char 2: 3-state ordered
  ordered3 <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3, 3, byrow = TRUE)

  # Char 1: A=0,B=0,C=1,D=1 -> 1
  # Char 2: A=0,B=0,C=2,D=2 -> 2
  res <- sankoff_test(
    tree4,
    n_states = c(2L, 3L),
    cost_matrices = list(fitch_cost(2), ordered3),
    tip_states = cbind(c(0, 0, 1, 1), c(0, 0, 2, 2))
  )
  expect_equal(res$score, 3)
  expect_equal(as.numeric(res$per_char), c(1, 2))
})


# ===== Larger tree ========================================================

test_that("5-tip tree scores correctly", {
  # (((t1,t2),t3),(t4,t5)), binary symmetric
  # t1=0, t2=0, t3=1, t4=1, t5=1 -> 1 step
  res <- sankoff_test(
    tree5,
    n_states = 2L,
    cost_matrices = list(fitch_cost(2)),
    tip_states = c(0, 0, 1, 1, 1)
  )
  expect_equal(res$score, 1)
})


# ===== X-transformation pattern (Goloboff et al. 2021) ====================

test_that("Asymmetric n+1:1 cost pattern matches expected score", {
  # 1 primary + 2 binary secondaries -> 5 states (absent + 4 present combos)
  # State 0 = absent
  # States 1-4 = (p00, p01, p10, p11) present combinations
  # gain = n+1 = 3, loss = 1
  # present -> present = Hamming distance of secondary states
  n <- 3  # gain cost
  cm <- matrix(0, 5, 5)
  for (i in 1:5) for (j in 1:5) {
    if (i == j) next
    if (i == 1) {
      cm[i, j] <- n  # absent -> any present = gain
    } else if (j == 1) {
      cm[i, j] <- 1  # any present -> absent = loss
    } else {
      # Hamming distance between binary encodings of (i-2) and (j-2)
      s1 <- c((i - 2) %/% 2, (i - 2) %% 2)
      s2 <- c((j - 2) %/% 2, (j - 2) %% 2)
      cm[i, j] <- sum(s1 != s2)
    }
  }

  # Tree ((A,B),(C,D)): A=absent(0), B=p00(1), C=p01(2), D=p11(4)
  # Optimal: root present, one gain; or root absent, multiple gains
  res <- sankoff_test(
    tree4,
    n_states = 5L,
    cost_matrices = list(cm),
    tip_states = c(0, 1, 2, 4),
    forced_root = 0L  # force root = absent (outgroup)
  )

  # Hand-computed: root=0(absent)
  # Need to gain at CD side and AB side independently
  # AB: A=absent(0), B=p00(1). Node AB: if present(p00), cost = 0(from B) + 1(loss from A->absent is wrong direction.. wait)
  # Actually: cost_matrix[parent_state][child_state]
  # If AB node = p00(1): min_t(cm[1][t]+cost[A][t]) = cm[1][0]+0 = 1 (present->absent=loss)
  #                       min_t(cm[1][t]+cost[B][t]) = cm[1][1]+0 = 0
  #                       AB cost for state 1: 1+0 = 1
  # If AB node = absent(0): min_t(cm[0][t]+cost[A][t]) = cm[0][0]+0 = 0
  #                          min_t(cm[0][t]+cost[B][t]) = cm[0][1]+0 = 3 (gain)
  #                          AB cost for state 0: 0+3 = 3
  # AB best = 1 (at state p00)

  # CD: C=p01(2), D=p11(4)
  # CD state p01(2): cost = cm[2][2]+0 + cm[2][4]+0 = 0 + 1 = 1
  # CD state p11(4): cost = cm[4][2]+0 + cm[4][4]+0 = 1 + 0 = 1
  # CD state absent(0): cost = cm[0][2]+0 + cm[0][4]+0 = 3 + 3 = 6
  # CD best = 1

  # Root forced at 0: cm[0][state_AB] + cost[AB][state_AB]
  #   state_AB=0: cm[0][0]+3 = 3
  #   state_AB=1: cm[0][1]+1 = 3+1 = 4
  #   best from AB: 3
  # cm[0][state_CD] + cost[CD][state_CD]
  #   state_CD=0: cm[0][0]+6 = 6
  #   state_CD=2: cm[0][2]+1 = 3+1 = 4
  #   state_CD=4: cm[0][4]+1 = 3+1 = 4
  #   best from CD: 4
  # Root cost = 3 + 4 = 7

  expect_equal(res$score, 7)
})
