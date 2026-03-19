#' Select a topologically diverse subset of trees
#'
#' Selects `n` trees from a `multiPhylo` object that are as
#' topologically distinct from one another as possible, using greedy
#' furthest-point (maximin) selection.  This is useful when a search
#' returns many most-parsimonious trees and downstream analyses
#' (consensus, tree-space visualization) need a manageable but
#' diverse subset.
#'
#' Uniform random subsampling of MPTs is misleading: the number of trees
#' in a region of tree space reflects the density of the parsimony
#' landscape, not the likelihood or support for that topology.
#' A random draw over-represents topologies that happen to sit on broad
#' plateaux and under-represents isolated optima.
#' `WideSample()` avoids this bias by selecting for topological spread
#' rather than frequency.
#'
#' @param trees A `multiPhylo` object.
#' @param n Integer: number of trees to retain. If `n >= length(trees)`, all
#'   trees are returned unchanged.
#' @param method Character: subsampling strategy.
#'   \describe{
#'     \item{`"maximin"`}{(Default.) Greedy furthest-point selection: start
#'       from a randomly chosen tree, then iteratively add the tree whose
#'       minimum distance to the already-selected set is largest. This
#'       maximizes topological spread.}
#'     \item{`"random"`}{Simple random sample without replacement.}
#'   }
#' @param distance Either:
#'   \itemize{
#'     \item A function that accepts a `multiPhylo` and returns a `dist`
#'       object (default: [TreeDist::ClusteringInfoDistance()]).
#'     \item A pre-computed `dist` object whose size matches `length(trees)`.
#'   }
#'   Ignored when `method = "random"`.
#'
#' @return A `multiPhylo` object of length `min(n, length(trees))`.
#'   Attributes of the input (e.g. `score`, `hits_to_best`) are preserved.
#'
#' @examples
#' library("TreeTools")
#' trees <- as.phylo(0:99, nTip = 8)
#' sub10 <- WideSample(trees, 10)
#' length(sub10)  # 10
#'
#' # Pre-computed distances
#' dists <- TreeDist::ClusteringInfoDistance(trees)
#' sub5 <- WideSample(trees, 5, distance = dists)
#'
#' @template MRS
#' @family tree scoring
#' @importFrom TreeDist ClusteringInfoDistance
#' @export
WideSample <- function(
    trees,
    n,
    method = c("maximin", "random"),
    distance = TreeDist::ClusteringInfoDistance
) {
  method <- match.arg(method)
  if (!inherits(trees, "multiPhylo")) {
    stop("`trees` must be a multiPhylo object")
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("`n` must be a single non-negative integer")
  }

  nTrees <- length(trees)
  if (n >= nTrees) {
    return(trees)
  }
  if (n == 0L) {
    return(.SubsetMultiPhylo(trees, integer(0)))
  }

  idx <- switch(method,
    random = sample.int(nTrees, n),
    maximin = .MaximinSubsample(trees, n, distance)
  )

  .SubsetMultiPhylo(trees, sort(idx))
}

#' @keywords internal
.MaximinSubsample <- function(trees, n, distance) {
  nTrees <- length(trees)
  if (inherits(distance, "dist")) {
    d <- as.matrix(distance)
    if (nrow(d) != nTrees) {
      stop(
        "`distance` has ", nrow(d), " entries but `trees` has ",
        nTrees, " trees"
      )
    }
  } else if (is.function(distance)) {
    d <- as.matrix(distance(trees))
  } else {
    stop("`distance` must be a function or a dist object")
  }

  selected <- integer(n)
  selected[1L] <- sample.int(nTrees, 1L)

  # min_dist[i] = min distance from tree i to any selected tree
  min_dist <- d[, selected[1L]]

  for (k in seq_len(n - 1L) + 1L) {
    # Zero out already-selected trees so they can't be picked again
    min_dist[selected[seq_len(k - 1L)]] <- -Inf
    selected[k] <- which.max(min_dist)
    min_dist <- pmin(min_dist, d[, selected[k]])
  }

  selected
}

#' Subset a multiPhylo preserving attributes
#' @keywords internal
.SubsetMultiPhylo <- function(trees, idx) {
  saved <- attributes(trees)
  result <- trees[idx]
  # Restore non-standard attributes (e.g. score, hits_to_best)
  standard <- c("names", "class")
  for (nm in setdiff(names(saved), standard)) {
    attr(result, nm) <- saved[[nm]]
  }
  result
}
