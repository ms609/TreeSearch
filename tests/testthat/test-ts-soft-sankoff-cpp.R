# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Compiled soft-Sankoff prototype (src/ts_soft_sankoff.{h,cpp}) against the
# verified pure-R reference and against two oracles that are independent of the
# soft path entirely.
#
# Plan: dev/plans/2026-08-01-soft-sankoff-temperature-dial.md
# R reference: tests/testthat/helper-soft-sankoff.R
#
# Matching the R reference is necessary but not sufficient: a shared conceptual
# error would pass such a test on both sides.  So the T = 1 case is also checked
# against FelsensteinLogLik() and phangorn::pml(), and the T -> 0 case against
# the pre-existing compiled hard kernel ts_sankoff_test().
#
# Every score comparison also asserts finiteness separately.  An out-of-range
# index or a missing per-branch matrix in this kernel yields Inf rather than an
# error, and expect_equal(Inf, Inf) passes.

RandomRootedCpp <- function(nTip) {
  TreeTools::Preorder(ape::rtree(nTip, rooted = TRUE, br = NULL))
}

RandomCostCpp <- function(nStates) {
  cost <- matrix(round(stats::runif(nStates^2, 0.5, 4), 2), nStates, nStates)
  diag(cost) <- 0
  cost
}

# Wrapper mirroring SoftSankoffScore()'s signature, so the two can be compared
# argument for argument.  `cost` is either a single matrix (shared by every
# edge) or a list indexed by child node, exactly as in the R reference.
CppSoftScore <- function(edge, tipCosts, cost, temperature, rootCost = NULL,
                         nRep = 1L) {
  nStates <- ncol(tipCosts)
  perBranch <- is.list(cost)
  storage.mode(edge) <- "integer"
  result <- TreeSearch:::ts_soft_sankoff_test(
    edge = edge,
    n_tip = nrow(tipCosts),
    tip_costs = list(tipCosts),
    # Ignored for a character that supplies branch_costs, but the binding still
    # validates its dimensions, so it must be conformable.
    cost_matrices = list(if (perBranch) matrix(0, nStates, nStates) else cost),
    temperature = temperature,
    root_costs = if (is.null(rootCost)) NULL else list(rootCost),
    branch_costs = if (perBranch) list(cost) else NULL,
    n_rep = as.integer(nRep)
  )
  result[["score"]]
}

# Assert two scores agree AND that neither is a silent infinity.
ExpectFiniteEqual <- function(observed, expected, tolerance = 1e-10) {
  expect_true(is.finite(observed))
  expect_true(is.finite(expected))
  expect_equal(observed, expected, tolerance = tolerance)
}

test_that("compiled soft score matches the R reference across T", {
  set.seed(101)
  for (rep in 1:12) {
    nTip <- sample(5:14, 1)
    nStates <- sample(2:5, 1)
    edge <- RandomRootedCpp(nTip)[["edge"]]
    states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
    tc <- TipCosts(states, nStates)
    cost <- RandomCostCpp(nStates)

    for (temperature in c(0, 1e-6, 0.02, 0.25, 1, 3)) {
      ExpectFiniteEqual(
        CppSoftScore(edge, tc, cost, temperature),
        SoftSankoffScore(edge, tc, cost, temperature)
      )
    }
  }
})

test_that("compiled kernel handles ambiguous tips (Inf costs) at small T", {
  # Infinite tip costs are the normal case, not an edge case: TipCosts() writes
  # 0 for observed and Inf for everything else.  A missing all-Inf guard in the
  # shifted log-sum-exp surfaces here as NaN, and a partially ambiguous tip is
  # the configuration that mixes finite and infinite entries in one reduction.
  set.seed(102)
  for (rep in 1:8) {
    nTip <- sample(6:12, 1)
    nStates <- 4
    edge <- RandomRootedCpp(nTip)[["edge"]]
    states <- replicate(nTip, sample(seq_len(nStates),
                                     sample(1:3, 1)), simplify = FALSE)
    tc <- TipCosts(states, nStates)
    cost <- RandomCostCpp(nStates)

    for (temperature in c(0, 0.02, 1)) {
      observed <- CppSoftScore(edge, tc, cost, temperature)
      expect_false(is.na(observed))
      ExpectFiniteEqual(observed,
                        SoftSankoffScore(edge, tc, cost, temperature))
    }
  }

  # A tip with every state impossible makes the whole score infinite, and must
  # do so without producing NaN from inf - inf.
  nStates <- 3
  edge <- RandomRootedCpp(6)[["edge"]]
  tc <- TipCosts(replicate(6, 1L, simplify = FALSE), nStates)
  tc[3, ] <- Inf
  cost <- RandomCostCpp(nStates)
  observed <- CppSoftScore(edge, tc, cost, 0.5)
  expect_false(is.na(observed))
  expect_identical(observed, Inf)
})

test_that("compiled T -> 0 reproduces the compiled hard Sankoff kernel", {
  # ts_sankoff_test() is an independent implementation: it does a hard min at
  # every node and at the root, with no temperature branch anywhere.
  #
  # It takes 0-based tip states (the guard in src/ts_rcpp.cpp is
  # `state >= 0 && state < ns_ch`); an out-of-range index there does not error,
  # it leaves every state at INF and returns Inf.  Hence the finiteness check.
  set.seed(103)
  for (rep in 1:8) {
    nTip <- sample(6:12, 1)
    nStates <- sample(2:4, 1)
    edge <- RandomRootedCpp(nTip)[["edge"]]
    states <- sample(seq_len(nStates), nTip, replace = TRUE)
    cost <- RandomCostCpp(nStates)

    hard <- TreeSearch:::ts_sankoff_test(
      edge = edge,
      n_states_r = as.integer(nStates),
      cost_matrices_r = list(cost),
      tip_states_r = matrix(as.integer(states) - 1L, ncol = 1L),
      forced_root_r = -1L
    )
    ExpectFiniteEqual(
      CppSoftScore(edge, TipCosts(as.list(states), nStates), cost, 0),
      as.numeric(hard[["score"]])
    )
  }
})

test_that("compiled T = 1 with per-branch -log(P) reproduces pruning", {
  # The load-bearing identity, checked against an oracle written in probability
  # space with no soft-min in it at all.  If this fails the prototype is wrong.
  set.seed(104)
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
    rootCost <- -log(rootFreq)

    observed <- CppSoftScore(edge, tc, cost, 1, rootCost = rootCost)
    ExpectFiniteEqual(observed,
                      -FelsensteinLogLik(edge, exp(-tc), P, rootFreq),
                      tolerance = 1e-12)
    ExpectFiniteEqual(observed,
                      SoftSankoffScore(edge, tc, cost, 1, rootCost = rootCost),
                      tolerance = 1e-12)
  }
})

test_that("compiled T = 1 matches phangorn::pml", {
  skip_if_not_installed("phangorn")
  set.seed(105)
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
  cost <- vector("list", max(edge))
  for (e in seq_len(nrow(edge))) {
    child <- edge[e, 2]
    bl <- rooted[["edge.length"]][e]
    cost[[child]] <- -log(MkTransition(2, if (is.na(bl)) 0 else bl))
  }
  lookup <- match(rooted[["tip.label"]], names(states))
  tc <- TipCosts(lapply(states[lookup],
                        function(s) if (s == "0") 1L else 2L), 2)

  observed <- -CppSoftScore(edge, tc, cost, 1, rootCost = -log(c(.5, .5)))
  ExpectFiniteEqual(observed, target, tolerance = 1e-8)
})

test_that("compiled score is monotone non-increasing in temperature", {
  set.seed(106)
  nTip <- 14
  nStates <- 3
  edge <- RandomRootedCpp(nTip)[["edge"]]
  states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
  tc <- TipCosts(states, nStates)
  cost <- matrix(1, nStates, nStates)
  diag(cost) <- 0

  grid <- c(0, 10^seq(-4, 0.5, length.out = 10))
  scores <- vapply(grid, function(tt) CppSoftScore(edge, tc, cost, tt),
                   numeric(1))
  expect_true(all(is.finite(scores)))
  expect_true(all(diff(scores) <= 1e-9))
})

test_that("compiled multi-character total is the sum of per-character scores", {
  set.seed(107)
  nTip <- 10
  edge <- RandomRootedCpp(nTip)[["edge"]]
  storage.mode(edge) <- "integer"
  nStatesVec <- c(2L, 3L, 5L, 2L)
  tipCostList <- lapply(nStatesVec, function(ns) {
    TipCosts(replicate(nTip, sample(seq_len(ns), 1), simplify = FALSE), ns)
  })
  costList <- lapply(nStatesVec, RandomCostCpp)

  result <- TreeSearch:::ts_soft_sankoff_test(
    edge = edge, n_tip = nTip,
    tip_costs = tipCostList, cost_matrices = costList,
    temperature = 0.3
  )
  expect_true(all(is.finite(result[["per_char"]])))
  expect_equal(result[["score"]], sum(result[["per_char"]]))

  # Characters of differing state counts share one strided tip-cost buffer
  # (stride = n_chars * max_states); a striding error would show up here as one
  # character reading another's costs.
  reference <- vapply(seq_along(nStatesVec), function(ch) {
    SoftSankoffScore(edge, tipCostList[[ch]], costList[[ch]], 0.3)
  }, numeric(1))
  expect_equal(result[["per_char"]], reference)
})

test_that("n_rep repeats the pass without changing the answer", {
  # n_rep exists so that Gate B times the kernel rather than R-side marshalling.
  # It must be a pure repetition: no accumulation, no state carried between
  # passes.
  set.seed(108)
  nTip <- 9
  nStates <- 3
  edge <- RandomRootedCpp(nTip)[["edge"]]
  tc <- TipCosts(replicate(nTip, sample(seq_len(nStates), 1),
                           simplify = FALSE), nStates)
  cost <- RandomCostCpp(nStates)

  once <- CppSoftScore(edge, tc, cost, 0.4, nRep = 1L)
  expect_true(is.finite(once))
  expect_identical(CppSoftScore(edge, tc, cost, 0.4, nRep = 25L), once)
})

test_that("the binding rejects trees it cannot score correctly", {
  nStates <- 2
  tc <- TipCosts(replicate(5, 1L, simplify = FALSE), nStates)
  cost <- matrix(c(0, 1, 1, 0), 2, 2)

  # A polytomy would be scored correctly by the R reference but not by a
  # binary-only kernel, so it must error rather than silently disagree.
  polytomy <- ape::read.tree(text = "((a,b,c),(d,e));")
  expect_error(CppSoftScore(polytomy[["edge"]], tc, cost, 0.5),
               "rooted binary tree")

  # An edge matrix with the right edge count but a three-child node: node 5 is
  # the root with children {6, 1, 2}, node 6 has the single child 7.
  threeChild <- cbind(c(5L, 5L, 5L, 6L, 7L, 7L),
                      c(6L, 1L, 2L, 7L, 3L, 4L))
  expect_error(
    CppSoftScore(threeChild, TipCosts(replicate(4, 1L, simplify = FALSE),
                                      nStates), cost, 0.5),
    "more than two children"
  )

  # Negative temperature is not a soft-min.
  edge <- RandomRootedCpp(5)[["edge"]]
  expect_error(CppSoftScore(edge, tc, cost, -1), "non-negative")
})
