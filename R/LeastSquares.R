# Least-squares distance tree fitting and search.
#
# A sibling to MaximizeParsimony() that optimises a least-squares fit to a
# target distance matrix instead of a parsimony score, using the same fast
# C++ rearrangement kernel (NNI + SPR).  Built for Lapointe & Cucumel's (1997)
# average consensus procedure, where the averaged patristic distance matrix is
# generally non-additive and the best-fitting topology must be found
# heuristically.

# Internal: coerce `dist` to a labelled symmetric matrix.
.LSMatrix <- function(dist) {
  D <- as.matrix(dist)
  if (nrow(D) != ncol(D)) {
    stop("`dist` must be a square distance matrix or a `dist` object")
  }
  labs <- rownames(D)
  if (is.null(labs)) labs <- colnames(D)
  if (is.null(labs)) {
    stop("`dist` must carry tip labels (row/column names)")
  }
  dimnames(D) <- list(labs, labs)
  D
}

# Internal: per-pair weight matrix, or NULL for unit weights.
# `weight`: NULL, "fm" (Fitch-Margoliash 1/D^2), or a numeric matrix.
.LSWeight <- function(weight, D) {
  if (is.null(weight)) return(NULL)
  if (is.character(weight)) {
    weight <- match.arg(weight, c("fm", "none"))
    if (weight == "none") return(NULL)
    W <- matrix(0, nrow(D), ncol(D), dimnames = dimnames(D))
    nz <- D != 0
    W[nz] <- 1 / (D[nz]^2)
    return(W)
  }
  W <- as.matrix(weight)
  if (!identical(dim(W), dim(D))) {
    stop("`weight` matrix must have the same dimensions as `dist`")
  }
  dimnames(W) <- dimnames(D)
  W
}

# Internal: prepare a starting tree for the C++ kernel — a rooted binary tree
# in canonical TreeSearch numbering (root = nTip + 1).  Returns the prepared
# tree; align the distance matrix to its `tip.label` order before calling C++.
.LSPrepTree <- function(tree, labs) {
  if (!inherits(tree, "phylo")) stop("starting tree must be a `phylo` object")
  tree <- TreeTools::KeepTip(tree, labs)
  # The kernel needs a *rooted* binary tree (n - 1 internal nodes, 2n - 2
  # edges).  Neighbour-joining and unrooted inputs are binary but unrooted, so
  # test for rootedness too; multi2di resolves any basal polytomy, rooting the
  # tree, and is a no-op on a tree that is already rooted and binary.
  if (!ape::is.binary(tree) || !ape::is.rooted(tree)) {
    tree <- ape::multi2di(tree, random = FALSE)
  }
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree[["tip.label"]])
  if (nrow(tree[["edge"]]) != 2L * nTip - 2L) {
    stop("Could not coerce starting tree to a rooted binary form")  # nocov
  }
  tree
}

# Internal: build the returned tree from a rooted binary edge matrix and fitted
# branch lengths.  Always constructs a *fresh* phylo (no inherited attributes)
# and strips any "order" attribute after unrooting: a stale order attribute
# (e.g. "preorder" from TreeTools::Preorder) makes ape's C routines, including
# cophenetic(), read the edge matrix in the wrong order and corrupt memory.
.LSFinalize <- function(edge, edgeLength, rss, tipLabels) {
  nTip <- length(tipLabels)
  out <- structure(
    list(edge = edge,
         edge.length = edgeLength,
         Nnode = nTip - 1L,
         tip.label = tipLabels),
    class = "phylo"
  )
  # phangorn convention: report the unrooted tree.  The kernel returns the two
  # root edges as (length, 0), so unrooting sums them to the true branch length.
  out <- ape::unroot(out)
  attr(out, "order") <- NULL
  attr(out, "RSS") <- rss
  out
}

#' Fit branch lengths to a distance matrix on a fixed topology
#'
#' Fits branch lengths on a fixed tree topology that minimise the (optionally
#' weighted) least-squares discrepancy between the tree's patristic distances
#' and a target distance matrix, using the package's C++ kernel.  This is the
#' fixed-topology counterpart of [`LeastSquaresTree()`], and the direct analogue
#' of [phangorn::nnls.tree()].
#'
#' @param tree A bifurcating tree of class \code{\link[ape]{phylo}}.  Edge
#' lengths, if any, are ignored and refitted.
#' @param dist A distance matrix (object of class \code{\link[stats]{dist}} or a
#' symmetric matrix with tip labels) over the tips of `tree`.
#' @param method Either `"nnls"` (non-negative least squares; branch lengths are
#' constrained to be \eqn{\ge 0}, matching [phangorn::nnls.tree()] and Lapointe
#' & Cucumel) or `"ols"` (ordinary least squares; faster, closed form, but may
#' return negative lengths).
#' @param weight Optional weighting of the residuals.  `NULL` (default) gives
#' unweighted least squares; `"fm"` applies Fitch-Margoliash weights
#' \eqn{1 / D_{ij}^2}; a numeric matrix supplies custom per-pair weights.
#'
#' @return The input `tree`, returned **unrooted**, with `edge.length` set to the
#' fitted branch lengths and an attribute `"RSS"` giving the residual sum of
#' squares.
#'
#' @examples
#' tree <- ape::rtree(8)
#' D <- cophenetic(tree)
#' fit <- LeastSquaresFit(tree, D)
#' attr(fit, "RSS")  # ~ 0: D is additive on this topology
#'
#' @seealso [`LeastSquaresTree()`] to search topologies; [phangorn::nnls.tree()].
#' @template MRS
#' @family least-squares functions
#' @importFrom stats cophenetic
#' @export
LeastSquaresFit <- function(tree, dist, method = c("nnls", "ols"),
                            weight = NULL) {
  method <- match.arg(method)
  D <- .LSMatrix(dist)
  W <- .LSWeight(weight, D)
  prepped <- .LSPrepTree(tree, rownames(D))
  labs <- prepped[["tip.label"]]
  Dord <- D[labs, labs, drop = FALSE]
  Word <- if (is.null(W)) NULL else W[labs, labs, drop = FALSE]
  methodCode <- if (method == "ols") 0L else 1L

  res <- ts_ls_fit(prepped[["edge"]], Dord, Word, methodCode)
  if (!isTRUE(res[["ok"]])) {
    warning("Least-squares solve was singular; results may be unreliable")
  }
  .LSFinalize(prepped[["edge"]], res[["edge_length"]], res[["rss"]], labs)
}

#' Find the least-squares-optimal tree for a distance matrix
#'
#' Searches tree topologies for the one whose patristic distances best fit a
#' target distance matrix under a least-squares criterion, fitting branch
#' lengths on each candidate and minimising the residual sum of squares.  The
#' heuristic uses the package's optimised C++ kernel, alternating \acronym{NNI}
#' and \acronym{SPR} rearrangements, exactly as the parsimony search does — but
#' driven by the least-squares score rather than tree length.
#'
#' This implements the topology-search step of Lapointe & Cucumel's (1997)
#' average consensus procedure, in which an averaged (and generally
#' non-additive) patristic distance matrix is fit by a Fitch-Margoliash
#' least-squares tree.
#'
#' @inheritParams LeastSquaresFit
#' @param dist A distance matrix (object of class \code{\link[stats]{dist}} or a
#' symmetric matrix with tip labels).
#' @param tree Optional starting point: a single \code{\link[ape]{phylo}} tree,
#' a list of trees (\code{multiPhylo}), or `NULL` (the default) to start from the
#' neighbour-joining tree of `dist`.  When several trees are supplied the search
#' is run from each and the best-fitting result is returned.
#' @param maxHits Integer; during hill-climbing, the number of equally-scoring
#' rearrangements to accept before moving on (helps traverse plateaux).
#' @param spr Logical; if `TRUE` (default) interleave \acronym{SPR} sweeps with
#' \acronym{NNI}, otherwise use \acronym{NNI} only (faster, more local).
#'
#' @return The best-fitting tree found, returned **unrooted**, with fitted
#' `edge.length` and an attribute `"RSS"` giving its residual sum of squares.
#'
#' @examples
#' set.seed(1)
#' trueTree <- ape::rtree(10)
#' D <- cophenetic(trueTree)        # additive: the generating tree fits exactly
#' found <- LeastSquaresTree(D)
#' attr(found, "RSS")               # ~ 0
#'
#' @seealso [`LeastSquaresFit()`] for fixed-topology fitting;
#' [`MaximizeParsimony()`] for the parsimony analogue.
#' @template MRS
#' @references \insertRef{LapointeCucumel1997}{TreeSearch}
#' @family least-squares functions
#' @importFrom stats cophenetic
#' @export
LeastSquaresTree <- function(dist, tree = NULL, method = c("nnls", "ols"),
                             weight = NULL, maxHits = 1L, spr = TRUE) {
  method <- match.arg(method)
  methodCode <- if (method == "ols") 0L else 1L
  D <- .LSMatrix(dist)
  W <- .LSWeight(weight, D)
  labs <- rownames(D)
  nTip <- length(labs)
  if (nTip < 4L) {
    stop("Least-squares tree search needs at least four tips")
  }

  starts <- if (is.null(tree)) {
    list(ape::nj(stats::as.dist(D)))
  } else if (inherits(tree, "phylo")) {
    list(tree)
  } else {
    # multiPhylo, possibly stored in compressed (.compressTipLabel) form where
    # components carry no `tip.label`.  Index with `[[`, whose multiPhylo method
    # restores the shared labels; `as.list()` would bypass it and yield
    # label-less trees.
    lapply(seq_along(tree), function(i) tree[[i]])
  }

  best <- NULL
  bestRSS <- Inf
  for (start in starts) {
    prepped <- .LSPrepTree(start, labs)
    tl <- prepped[["tip.label"]]
    Dord <- D[tl, tl, drop = FALSE]
    Word <- if (is.null(W)) NULL else W[tl, tl, drop = FALSE]

    res <- ts_ls_search(prepped[["edge"]], Dord, Word, methodCode,
                        as.integer(maxHits), isTRUE(spr))
    # Keep the first result unconditionally so a singular fit (RSS = Inf, e.g.
    # a weighting that leaves a branch unidentifiable) still yields a tree
    # rather than NULL; better fits replace it.
    if (is.null(best) || res[["rss"]] < bestRSS) {
      bestRSS <- res[["rss"]]
      best <- .LSFinalize(res[["edge"]], res[["edge_length"]], res[["rss"]], tl)
    }
  }
  if (!is.finite(bestRSS)) {
    warning("Least-squares fit was singular for every starting tree; ",
            "branch lengths are unreliable. Check for zero weights/distances.")
  }
  best
}
