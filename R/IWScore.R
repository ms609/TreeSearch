#' @template pointlessDots
#' @rdname TreeLength
#' @export
IWScore <- function (tree, dataset, concavity = 10L, ...) {
  .Deprecated('TreeLength')
  TreeLength(tree, dataset, concavity)
}

#' @describeIn TreeSearch Search using implied weights.
#' @template concavityParam
#' @export
IWTreeSearch <- function (tree, dataset, concavity = 10L, 
                          EdgeSwapper = RootedTBR,
                          maxIter = 100L, maxHits = 20L,
                          verbosity = 1L, ...) {
  .Deprecated("MaximizeParsimony") # Retained as template, for now.
  #TODO move all these functions to a vignette.
  if (!inherits(dataset, 'phyDat')) {
    stop("Unrecognized dataset class; should be phyDat, not ",
         class(dataset), '.')
  }
  if (!('min.length' %in% names(attributes(dataset)))) {
    dataset <- PrepareDataIW(dataset)
  }
  at <- attributes(dataset)
  
  TreeSearch(tree, dataset, nChar=at$nr, weight=at$weight,
             minLength=at$min.length, concavity = concavity,
             InitializeData = IWInitMorphy,
             CleanUpData = IWDestroyMorphy,
             TreeScorer = IWScoreMorphy,
             EdgeSwapper = EdgeSwapper,
             maxIter = maxIter, maxHits = maxHits, verbosity = verbosity, ...)
}
