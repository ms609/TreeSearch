#' Rearrange edges of a phylogenetic tree
#' 
#' `RearrangeEdges()` performs the specified edge rearrangement on a matrix 
#' that corresponds to the edges of a phylogenetic tree, returning the score of
#' the new tree.
#' Will generally be called from within a tree search function.
#' 
#' @details `RearrangeTree()` performs one tree rearrangement of a
#'  specified type, and returns the score of the tree (with the given dataset).
#'  It also reports the number of times that this score was hit in the 
#'  current function call.
#' 
#' @template treeParent
#' @template treeChild
#' @param dataset Third argument to pass to \code{TreeScorer}.
#' @template treeScorerParam
#' @param scoreToBeat Double giving score of input tree.
#' @param hits Integer giving number of times the input tree has already been hit.
#' @template EdgeSwapperParam
## @param  minScore trees longer than \code{minScore}, probably the score of the best previously known tree,
##     will be discarded;
## @param returnSingle returns all trees if `FALSE` or a randomly selected tree if `TRUE`.
#' @param iter iteration number of calling function, for reporting to user only.
#' @template verbosityParam
#' @template treeScorerDots
#'
#' @author Martin R. Smith
#'
#' @template returnEdgeList
#' 
#' @examples
#' data('Lobo', package='TreeTools')
#' tree <- TreeTools::NJTree(Lobo.phy)
#' edge <- tree$edge
#' parent <- edge[, 1]
#' child <- edge[, 2]
#' dataset <- PhyDat2Morphy(Lobo.phy)
#' RearrangeEdges(parent, child, dataset, EdgeSwapper = RootedNNISwap)
#' # Remember to free memory:
#' dataset <- UnloadMorphy(dataset)
#' @export
RearrangeEdges <- function (parent, child, dataset, TreeScorer = MorphyLength,
                            EdgeSwapper, 
                            scoreToBeat = TreeScorer(parent, child, dataset, ...),
                            iter = '?', hits = 0L, verbosity = 0L, ...) {
  eps <- .Machine$double.eps ^ 0.5
  rearrangedEdges <- EdgeSwapper(parent, child)
  if (is.list(rearrangedEdges[[1]])) {
    # Then we've been sent a list of possible trees
    candidateScores <- vapply(rearrangedEdges, function (edges) {
      TreeScorer(edges[[1]], edges[[2]], dataset, ...)
    } , double(1))
    candidateScore <- min(candidateScores)
    best <- candidateScores == candidateScore
    nBest <- sum(best)
    if (candidateScore > (scoreToBeat + eps)) {
      if (verbosity > 3L) {                                                     # nocov start
        message("    . Iteration ", iter, 
                                  ' - Rearranged tree score ', candidateScore, 
                              " > target ", scoreToBeat)
      }                                                                         # nocov end
    } else if (candidateScore + eps > scoreToBeat) { # i.e. scores are equal
      hits <- hits + nBest
      if (verbosity > 2L) {                                                     # nocov start
        message("    - Iteration ", iter, " - Best score ", scoreToBeat, 
                " found again ", nBest, " times; now found ", hits, " times.")
      }                                                                         # nocov end
    } else {
      hits <- nBest
      if (verbosity > 1L) {                                                     # nocov start
        message("    * Iteration ", iter, " - New best score ", candidateScore,
                " found on ", hits, " trees.")
      }                                                                         # nocov end
    }
    rearrangedEdges <- rearrangedEdges[[SampleOne(which(best), nBest)]]
  } else {
    candidateScore <- TreeScorer(rearrangedEdges[[1]], rearrangedEdges[[2]], dataset, ...)
    if (candidateScore > (scoreToBeat + eps)) {
      if (verbosity > 3L) {                                                     # nocov start
        message("    . Iteration ", iter, ' - Rearranged tree score ',
                signif(candidateScore, 6), " > target ", signif(scoreToBeat, 6))
      }                                                                         # nocov end
    } else if (candidateScore + eps > scoreToBeat) { # i.e. scores are equal
      hits <- hits + 1L
      if (verbosity > 2L) {                                                     # nocov start
        message("    - Iteration ", iter, " - Best score ", 
                signif(scoreToBeat, 6), " hit ", hits, " times.")
      }                                                                         # nocov end
    } else {
      hits <- 1L
      if (verbosity > 1L) {                                                     # nocov start
        message("    * Iteration ", iter, " - New best score ",
                signif(candidateScore, 6), " found on ", hits, " trees.")
      }                                                                         # nocov end
    }
  }
  rearrangedEdges[3:4] <- c(candidateScore, hits)
  # Return:
  rearrangedEdges
}

#' Check that all nodes in a tree are bifurcating.
#' 
#' @template treeParent
#' 
#' @return Returns `NULL`, but will `stop` with an error message if a tree
#' does not appear to be bifurcating.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
StopUnlessBifurcating <- function (parent) {
  if (!all(table(parent) == 2L)) stop ("Tree must be bifurcating; try collapse.singles or multi2di.")
}
