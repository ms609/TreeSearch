# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Exploratory tests for the soft-Sankoff temperature dial.
# Plan: dev/plans/2026-08-01-soft-sankoff-temperature-dial.md
# Reference implementation: tests/testthat/helper-soft-sankoff.R
#
# These pin the two claims the whole idea rests on:
#   T -> 0  ==  hard weighted parsimony (and, with equal costs, Fitch)
#   T  = 1  ==  Felsenstein pruning, with cost = -log(P)
# If the second fails, the plan document is wrong and should be withdrawn
# rather than debugged.

RandomRooted <- function(nTip) {
  TreeTools::Preorder(ape::rtree(nTip, rooted = TRUE, br = NULL))
}

RandomCost <- function(nStates) {
  cost <- matrix(round(stats::runif(nStates^2, 0.5, 4), 2), nStates, nStates)
  diag(cost) <- 0
  cost
}

test_that("SoftMin() converges to min at the rate set by multiplicity", {
  x <- c(3.2, 1.7, 1.7, 9)
  expect_equal(SoftMin(x, 0), min(x))
  expect_true(all(vapply(10^seq(0, -6), function(tt) {
    SoftMin(x, tt) <= min(x) + 1e-12
  }, logical(1))))

  # softmin -> min - T * log(#minima).  That gap is the entropy of the
  # near-optimal set: the extent to which the criterion integrates over
  # reconstructions rather than optimising them.
  expect_equal((min(x) - SoftMin(x, 1e-7)) / 1e-7, log(2), tolerance = 1e-6)
  expect_equal((min(c(1, 1, 1, 5)) - SoftMin(c(1, 1, 1, 5), 1e-7)) / 1e-7,
               log(3), tolerance = 1e-6)

  expect_equal(SoftMin(c(Inf, Inf), 1), Inf)
  expect_equal(SoftMin(c(Inf, 2, Inf), 1e-6), 2)
})

test_that("T -> 0 reproduces hard Sankoff", {
  set.seed(11)
  for (rep in 1:10) {
    nTip <- sample(5:12, 1)
    nStates <- sample(2:5, 1)
    edge <- RandomRooted(nTip)[["edge"]]
    states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
    tc <- TipCosts(states, nStates)
    cost <- RandomCost(nStates)
    expect_equal(SoftSankoffScore(edge, tc, cost, temperature = 1e-8),
                 HardSankoffScore(edge, tc, cost),
                 tolerance = 1e-5)
  }
})

test_that("T -> 0 reproduces the compiled Sankoff kernel", {
  set.seed(12)
  for (rep in 1:6) {
    nTip <- sample(6:12, 1)
    nStates <- sample(2:4, 1)
    tree <- RandomRooted(nTip)
    edge <- tree[["edge"]]
    states <- sample(seq_len(nStates), nTip, replace = TRUE)
    cost <- RandomCost(nStates)

    # `ts_sankoff_test()` takes 0-based tip states (src/ts_rcpp.cpp: the guard
    # is `state >= 0 && state < ns_ch`).  The "1-based" in
    # R/recode_hierarchy.R:188 means present states begin at index 1 of a
    # 0-based array, state 0 being "absent"; an out-of-range index silently
    # leaves every state at INF rather than erroring.
    compiled <- TreeSearch:::ts_sankoff_test(
      edge = edge,
      n_states_r = as.integer(nStates),
      cost_matrices_r = list(cost),
      tip_states_r = matrix(as.integer(states) - 1L, ncol = 1L),
      forced_root_r = -1L
    )
    reference <- SoftSankoffScore(
      edge, TipCosts(as.list(states), nStates), cost, temperature = 1e-8
    )
    expect_equal(as.numeric(compiled[["score"]]), reference, tolerance = 1e-5)
  }
})

test_that("T = 1 with cost = -log(P) reproduces Felsenstein pruning", {
  set.seed(13)
  for (rep in 1:10) {
    nTip <- sample(5:12, 1)
    nStates <- sample(2:4, 1)
    tree <- TreeTools::Preorder(ape::rtree(nTip, rooted = TRUE))
    edge <- tree[["edge"]]

    P <- cost <- vector("list", max(edge))
    for (e in seq_len(nrow(edge))) {
      child <- edge[e, 2]
      P[[child]] <- MkTransition(nStates, tree[["edge.length"]][e])
      cost[[child]] <- -log(P[[child]])
    }
    states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
    tc <- TipCosts(states, nStates)
    rootFreq <- rep(1 / nStates, nStates)

    expect_equal(
      SoftSankoffScore(edge, tc, cost, temperature = 1,
                       rootCost = -log(rootFreq)),
      -FelsensteinLogLik(edge, exp(-tc), P, rootFreq),
      tolerance = 1e-12
    )
  }
})

test_that("T = 1 matches phangorn::pml on a two-state dataset", {
  skip_if_not_installed("phangorn")
  set.seed(14)
  nTip <- 8
  tree <- ape::rtree(nTip, rooted = FALSE)
  tree[["edge.length"]] <- round(stats::runif(nrow(tree[["edge"]]), .02, .4), 3)
  states <- sample(c("0", "1"), nTip, replace = TRUE)
  names(states) <- tree[["tip.label"]]
  dat <- phangorn::phyDat(as.matrix(states), type = "USER",
                          levels = c("0", "1"))
  target <- as.numeric(phangorn::pml(tree, dat)[["logLik"]])

  rooted <- ape::multi2di(ape::root(tree, outgroup = tree[["tip.label"]][1],
                                    resolve.root = TRUE))
  edge <- rooted[["edge"]]
  P <- cost <- vector("list", max(edge))
  for (e in seq_len(nrow(edge))) {
    child <- edge[e, 2]
    bl <- rooted[["edge.length"]][e]
    P[[child]] <- MkTransition(2, if (is.na(bl)) 0 else bl)
    cost[[child]] <- -log(P[[child]])
  }
  lookup <- match(rooted[["tip.label"]], names(states))
  tc <- TipCosts(lapply(states[lookup],
                        function(s) if (s == "0") 1L else 2L), 2)

  expect_equal(
    -SoftSankoffScore(edge, tc, cost, temperature = 1,
                      rootCost = -log(c(.5, .5))),
    target, tolerance = 1e-8
  )
})

test_that("equal costs at T -> 0 reproduce the Fitch score", {
  set.seed(15)
  for (rep in 1:8) {
    nTip <- sample(6:14, 1)
    nStates <- sample(2:4, 1)
    tree <- RandomRooted(nTip)
    states <- sample(seq_len(nStates) - 1L, nTip, replace = TRUE)
    cost <- matrix(1, nStates, nStates)
    diag(cost) <- 0

    soft <- SoftSankoffScore(tree[["edge"]],
                             TipCosts(as.list(states + 1L), nStates),
                             cost, temperature = 1e-8)
    expect_equal(soft, FitchScore(tree[["edge"]], as.list(states + 1L),
                                  nStates),
                 tolerance = 1e-5)

    # ... and the package's own Fitch, which is rooting-invariant.
    mat <- matrix(as.character(states), nrow = nTip,
                  dimnames = list(tree[["tip.label"]], NULL))
    dataset <- phangorn::phyDat(mat, type = "USER",
                                levels = as.character(seq_len(nStates) - 1L))
    expect_equal(soft, as.numeric(TreeLength(tree, dataset, concavity = Inf)),
                 tolerance = 1e-5)
  }
})

test_that("score is monotone non-increasing in temperature", {
  set.seed(16)
  nTip <- 12
  nStates <- 3
  edge <- RandomRooted(nTip)[["edge"]]
  states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
  tc <- TipCosts(states, nStates)
  cost <- matrix(1, nStates, nStates)
  diag(cost) <- 0

  grid <- c(0, 10^seq(-4, 0.5, length.out = 10))
  scores <- vapply(grid, function(tt) {
    SoftSankoffScore(edge, tc, cost, temperature = tt)
  }, numeric(1))
  expect_true(all(diff(scores) <= 1e-9))
  expect_equal(scores[1], SoftSankoffScore(edge, tc, cost, 0))
})

test_that("up-pass marginals are proper and respect observed tip states", {
  set.seed(17)
  nTip <- 10
  nStates <- 3
  edge <- RandomRooted(nTip)[["edge"]]
  states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
  tc <- TipCosts(states, nStates)
  cost <- matrix(1, nStates, nStates)
  diag(cost) <- 0

  marginals <- SoftSankoffMarginals(edge, tc, cost, temperature = 0.5)
  expect_equal(rowSums(marginals), rep(1, nrow(marginals)))

  # A tip's marginal must be an indicator of its observed state.
  for (tip in seq_len(nTip)) {
    expect_equal(marginals[tip, states[[tip]]], 1)
  }

  # Higher temperature spreads the internal-node marginals: mean entropy must
  # increase with T.  This is the property Step 3b would surface to users.
  Entropy <- function(temperature) {
    m <- SoftSankoffMarginals(edge, tc, cost, temperature)
    internal <- m[(nTip + 1):nrow(m), , drop = FALSE]
    mean(apply(internal, 1, function(p) {
      p <- p[p > 0]
      -sum(p * log(p))
    }))
  }
  expect_lt(Entropy(0.1), Entropy(1))
})
