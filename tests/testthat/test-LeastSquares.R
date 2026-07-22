library("TreeTools", quietly = TRUE)

# A reproducible additive matrix: patristic distances of a known tree.
additiveCase <- function(n, seed = 1) {
  set.seed(seed)
  truth <- ape::rtree(n, br = function(k) runif(k, 0.5, 3))
  list(truth = truth, D = cophenetic(truth))
}

# A reproducible non-additive matrix over n labelled tips.
nonAdditiveCase <- function(n, seed = 1) {
  set.seed(seed)
  labs <- paste0("t", seq_len(n))
  M <- matrix(0, n, n, dimnames = list(labs, labs))
  M[lower.tri(M)] <- runif(n * (n - 1) / 2, 1, 5)
  M + t(M)
}

sameTopology <- function(a, b) {
  isTRUE(TreeDist::RobinsonFoulds(a, b) == 0)
}

test_that("LeastSquaresFit fits additive matrices exactly", {
  for (n in c(5, 8, 11)) {
    case <- additiveCase(n)
    for (method in c("nnls", "ols")) {
      fit <- LeastSquaresFit(case$truth, case$D, method = method)
      expect_lt(attr(fit, "RSS"), 1e-8)
      # Fitted patristic distances reproduce the target.
      coph <- cophenetic(fit)[case$truth$tip.label, case$truth$tip.label]
      expect_equal(coph, case$D[case$truth$tip.label, case$truth$tip.label],
                   tolerance = 1e-6)
    }
  }
})

test_that("LeastSquaresFit matches phangorn::nnls.tree (RSS + branch lengths)", {
  skip_if_not_installed("phangorn")
  for (n in c(6, 9, 12)) {
    # Additive and non-additive targets on a fixed topology.
    set.seed(n)
    tree <- ape::rtree(n, br = function(k) runif(k, 0.2, 2))
    labs <- tree$tip.label
    for (case in c("additive", "nonadditive")) {
      D <- if (case == "additive") {
        cophenetic(tree)[labs, labs]
      } else {
        nonAdditiveCase(n, seed = n + 100)[labs, labs]
      }

      ph <- phangorn::nnls.tree(D, tree, method = "unrooted")
      rssPh <- attr(ph, "RSS")
      if (is.null(rssPh)) {                       # quadprog branch omits RSS
        cph <- cophenetic(ph)[labs, labs]
        rssPh <- sum((D[lower.tri(D)] - cph[lower.tri(cph)])^2)
      }

      fit <- LeastSquaresFit(tree, D, method = "nnls")
      # RSS agreement
      expect_equal(attr(fit, "RSS"), rssPh, tolerance = 1e-5)
      # Fitted distance agreement (parametrisation-independent branch lengths)
      cphMe <- cophenetic(fit)[labs, labs]
      cphPh <- cophenetic(ph)[labs, labs]
      expect_equal(cphMe[lower.tri(cphMe)], cphPh[lower.tri(cphPh)],
                   tolerance = 1e-5)
    }
  }
})

test_that("LeastSquaresFit OLS matches a direct normal-equation solve", {
  skip_if_not_installed("phangorn")
  n <- 10
  set.seed(3)
  tree <- ape::rtree(n)
  D <- nonAdditiveCase(n, seed = 5)[tree$tip.label, tree$tip.label]

  ut <- ape::unroot(tree)
  X <- as.matrix(phangorn::designTree(ut))      # unrooted path design matrix
  dm <- D[ut$tip.label, ut$tip.label]
  y <- dm[lower.tri(dm)]
  beta <- solve(crossprod(X), crossprod(X, y))
  rssDirect <- sum((y - X %*% beta)^2)

  fit <- LeastSquaresFit(tree, D, method = "ols")
  expect_equal(attr(fit, "RSS"), rssDirect, tolerance = 1e-5)
})

test_that("LeastSquaresTree recovers the generating tree from an additive matrix", {
  for (n in c(6, 8, 10, 12)) {
    case <- additiveCase(n, seed = n)
    RNGkind("Mersenne-Twister")

    # Neighbour-joining seed (exact on additive data): must hold the optimum.
    set.seed(1)
    fromNJ <- LeastSquaresTree(case$D)
    expect_lt(attr(fromNJ, "RSS"), 1e-6)
    expect_true(sameTopology(fromNJ, case$truth))

    # A deliberately poor random start: the search must climb to the optimum.
    set.seed(2)
    badStart <- ape::rtree(n, tip.label = case$truth$tip.label,
                           br = function(k) rep(1, k))
    fromBad <- LeastSquaresTree(case$D, tree = badStart)
    expect_lt(attr(fromBad, "RSS"), 1e-6)
    expect_true(sameTopology(fromBad, case$truth))
  }
})

test_that("LeastSquaresTree reduces RSS on a non-additive matrix", {
  D <- nonAdditiveCase(12, seed = 99)
  start <- ape::rtree(12, tip.label = rownames(D), br = function(k) rep(1, k))
  startRSS <- attr(LeastSquaresFit(start, D), "RSS")

  RNGkind("Mersenne-Twister")
  set.seed(1)
  found <- LeastSquaresTree(D, tree = start)
  expect_s3_class(found, "phylo")
  expect_false(ape::is.rooted(found))            # phangorn convention
  expect_lt(attr(found, "RSS"), startRSS)         # genuinely improved
  expect_equal(length(found$edge.length), nrow(found$edge))
})

test_that("OLS and NNLS agree when the OLS fit is already non-negative", {
  case <- additiveCase(9, seed = 4)
  ols <- LeastSquaresFit(case$truth, case$D, method = "ols")
  nnls <- LeastSquaresFit(case$truth, case$D, method = "nnls")
  cOls <- cophenetic(ols)[case$truth$tip.label, case$truth$tip.label]
  cNnls <- cophenetic(nnls)[case$truth$tip.label, case$truth$tip.label]
  expect_equal(cOls, cNnls, tolerance = 1e-6)
})

test_that("Fitch-Margoliash weighting matches a direct weighted OLS solve", {
  skip_if_not_installed("phangorn")
  n <- 9
  set.seed(7)
  tree <- ape::rtree(n)
  D <- nonAdditiveCase(n, seed = 21)[tree$tip.label, tree$tip.label]

  # Direct weighted OLS via the unrooted path design matrix (label-aligned,
  # so no dependence on phangorn's internal pair ordering).
  ut <- ape::unroot(tree)
  X <- as.matrix(phangorn::designTree(ut))
  labs <- ut$tip.label
  dm <- D[labs, labs]
  wm <- 1 / (dm^2); diag(wm) <- 0
  y <- dm[lower.tri(dm)]
  wv <- wm[lower.tri(wm)]
  Xw <- X * sqrt(wv)
  yw <- y * sqrt(wv)
  beta <- solve(crossprod(Xw), crossprod(Xw, yw))
  rssDirect <- sum((yw - Xw %*% beta)^2)

  fitW <- LeastSquaresFit(tree, D, method = "ols", weight = "fm")
  expect_equal(attr(fitW, "RSS"), rssDirect, tolerance = 1e-4)

  # Weighting genuinely changes the fitted distances relative to unweighted.
  fitU <- LeastSquaresFit(tree, D, method = "ols")
  cW <- cophenetic(fitW)[tree$tip.label, tree$tip.label]
  cU <- cophenetic(fitU)[tree$tip.label, tree$tip.label]
  expect_gt(max(abs(cW - cU)), 1e-6)
})

test_that("LeastSquaresTree accepts multiple starting trees (incl. compressed)", {
  case <- additiveCase(9, seed = 5)
  labs <- case$truth$tip.label
  set.seed(6)
  starts <- structure(list(
    ape::rtree(9, tip.label = labs, br = function(k) rep(1, k)),
    ape::rtree(9, tip.label = labs, br = function(k) rep(1, k))
  ), class = "multiPhylo")

  # Plain (uncompressed) multiPhylo: runs from each start, returns the best.
  set.seed(1)
  fromMulti <- LeastSquaresTree(case$D, tree = starts)
  expect_s3_class(fromMulti, "phylo")
  expect_lt(attr(fromMulti, "RSS"), 1e-6)
  expect_true(sameTopology(fromMulti, case$truth))

  # Compressed (.compressTipLabel) form: components carry no tip.label, so the
  # search must restore them via `[[`.  Result must be no worse.
  compressed <- ape::.compressTipLabel(starts)
  set.seed(1)
  fromCompressed <- LeastSquaresTree(case$D, tree = compressed)
  expect_s3_class(fromCompressed, "phylo")
  expect_lt(attr(fromCompressed, "RSS"), 1e-6)
  expect_true(sameTopology(fromCompressed, case$truth))
})

test_that("A rank-deficient (zero-weight) fit fails gracefully, not fatally", {
  # Zero weights / zero distances can leave a branch unidentifiable, making the
  # OLS normal equations singular.  The fit must signal failure and return a
  # finite, fully-sized length vector rather than reading off the end of an
  # empty branch_length buffer (which previously segfaulted).
  set.seed(2)
  tree <- ape::rtree(6)
  D <- cophenetic(tree)
  W <- matrix(1, 6, 6, dimnames = dimnames(D)); diag(W) <- 0
  W[1, 2] <- 0; W[2, 1] <- 0

  prepped <- .LSPrepTree(tree, rownames(D))
  tl <- prepped[["tip.label"]]
  res <- ts_ls_fit(prepped[["edge"]], D[tl, tl], W[tl, tl], 0L)  # 0L = OLS
  expect_false(isTRUE(res[["ok"]]))
  # Pin the failure contract exactly: a fully-sized, all-zero length vector and
  # RSS = Inf.  This detects a *non-crashing* reversion, which a soft
  # all(is.finite()) check could miss if an out-of-bounds read returned garbage.
  expect_equal(res[["edge_length"]], rep(0, nrow(prepped[["edge"]])))
  expect_equal(res[["rss"]], Inf)

  # User-facing path: FM weighting with a zero distance -> zero weight.  Must
  # warn and return a valid phylo, never crash.
  set.seed(3)
  tree2 <- ape::rtree(6)
  D2 <- cophenetic(tree2); D2[1, 2] <- 0; D2[2, 1] <- 0
  expect_warning(fit <- LeastSquaresFit(tree2, D2, method = "ols", weight = "fm"),
                 "singular")
  expect_s3_class(fit, "phylo")
})

test_that("LeastSquaresTree accepts dist objects and validates input", {
  case <- additiveCase(7, seed = 8)
  d <- stats::as.dist(case$D)
  set.seed(1)
  found <- LeastSquaresTree(d)
  expect_s3_class(found, "phylo")
  expect_lt(attr(found, "RSS"), 1e-6)

  # Unlabelled matrix is rejected.
  bare <- unname(case$D)
  expect_error(LeastSquaresTree(bare), "tip label")

  # Too few tips.
  tiny <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3,
                 dimnames = list(letters[1:3], letters[1:3]))
  expect_error(LeastSquaresTree(tiny), "four tips")
})
