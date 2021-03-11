#' @template pointlessDots
#' @rdname TreeLength
#' @export
IWScore <- function (tree, dataset, concavity = 10, ...) {
  .Deprecated('TreeLength')
  TreeLength(tree, dataset, concavity)
}

#' @describeIn ProfileScore Scorer for Implied Weighting dataset.
#' @template concavityParam
#' @param minLength Integer vector specifying the minimum length
#'                  possible for each character in `dataset`, perhaps calculated
#'                  using \code{\link{MinimumLength}}.
#'
#' @export
IWScoreMorphy <- function (parent, child, dataset, concavity = 10L, 
                           minLength = attr(dataset, 'min.length'), ...) {
  steps <- vapply(attr(dataset, 'morphyObjs'), MorphyLength,
                  parent = parent, child = child, integer(1))
  homoplasies <- steps - minLength
  fit <- homoplasies / (homoplasies + concavity)
  # Return:
  sum(fit * attr(dataset, 'weight'))
}

#' @rdname ProfileScore
#' @export
IWInitMorphy <- function (dataset) {
  attr(dataset, 'morphyObjs') <- 
    lapply(PhyToString(dataset, byTaxon = FALSE, useIndex = FALSE, 
                       concatenate = FALSE), 
           SingleCharMorphy)
  
  # Return:
  dataset
}


#' @describeIn TreeSearch Search using profile parsimony
#' @template concavityParam
#' @export
IWTreeSearch <- function (tree, dataset, concavity = 10, 
                          EdgeSwapper = RootedTBR,
                          maxIter = 100, maxHits = 20,
                          verbosity = 1, ...) {
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