#' @inheritParams TreeSearch
#' @inheritParams EdgeListSearch
#' @param dataset A `ParsimonyData` object from [`PrepareData()`].
#' @param maxIter Numeric specifying maximum number of iterations to perform in
#' tree search.
#' @param maxHits Numeric specifying maximum number of hits to accomplish in
#' tree search.
#' @param stopAtScore stop search as soon as this score is hit or beaten.
#' @param \dots further parameters to send to `TreeScorer()`
#'
#' @return `BootstrapTree()` returns a tree that is optimal under a random
#' sampling of the original characters.
#'
#' @rdname Ratchet
#' @export
BootstrapTree <- function (edgeList, dataset, EdgeSwapper = NNISwap,
                           maxIter, maxHits, verbosity = 1L,
                           stopAtPeak = FALSE, stopAtPlateau = 0L, ...) {
  startWeights <- dataset[["original_weight"]]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep.int(eachChar, startWeights)
  resampling <- tabulate(sample(deindexedChars, replace = TRUE),
                         length(startWeights))
  # R copy-on-modify: the caller's `dataset` is unchanged.
  dataset[["weight"]] <- as.integer(resampling)

  res <- EdgeListSearch(edgeList[1:2], dataset, EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        stopAtPeak = stopAtPeak, stopAtPlateau = stopAtPlateau,
                        verbosity = verbosity - 1L, ...)

  # Return:
  res[1:2]
}
