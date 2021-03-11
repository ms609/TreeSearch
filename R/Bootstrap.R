#' Ratchet bootstrapper
#' 
#' @template edgeListParam
#' @template morphyObjParam
#' @template EdgeSwapperParam
#' @param maxIter maximum number of iterations to perform in tree search
#' @param maxHits maximum number of hits to accomplish in tree search
#' @template stopAtPeakParam
#' @template stopAtPlateauParam
#' @template verbosityParam
#' @param \dots further parameters to send to `TreeScorer()`
#'
#' @return A tree that is optimal under a random sampling of the original characters
#' @export
MorphyBootstrap <- function (edgeList, morphyObj, EdgeSwapper = NNISwap, 
                             maxIter, maxHits, verbosity = 1L, 
                             stopAtPeak = FALSE, stopAtPlateau=0L, ...) {
  startWeights <- MorphyWeights(morphyObj)['exact', ]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep(eachChar, startWeights)
  resampling <- tabulate(sample(deindexedChars, replace = TRUE),
                         length(startWeights))
  errors <- vapply(eachChar, function (i) 
            mpl_set_charac_weight(i, resampling[i], morphyObj), integer(1))
  
  if (any(errors)) {
    stop ("Error resampling morphy object: ",
          mpl_translate_error(unique(errors[errors < 0L])))
  }
  if (mpl_apply_tipdata(morphyObj) -> error) {
    stop("Error applying tip data: ", mpl_translate_error(error))
  }
  
  res <- EdgeListSearch(edgeList[1:2], morphyObj, EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        stopAtPeak = stopAtPeak, stopAtPlateau = stopAtPlateau,
                        verbosity = verbosity - 1L, ...)
  errors <- vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
  if (any(errors)) {
    stop ("Error resampling morphy object: ",
          mpl_translate_error(unique(errors[errors < 0L])))
  }
  if (mpl_apply_tipdata(morphyObj) -> error) {
    stop("Error applying tip data: ", mpl_translate_error(error))
  }
  
  # Return:
  res[1:2]
}
