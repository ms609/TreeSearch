#' Implied weights parsimony Score
#'
#' Calculate a tree's Parsimony score with a given dataset using implied weights
#' (Goloboff 1997).
#'
#' @template treeParam
#' @param dataset Dataset of class \code{phyDat}.  The dataset should have been
#' prepared using \code{dataset <- \link{PrepareDataIW}(dataset)}; if this step
#' has not been completed, the dataset will be (time-consumingly) prepared
#' within the function.
#' In subsidiary functions, the dataset will have been initialized using 
#' \code{IWInitMorphy}, and must be destroyed using \code{IWDestroyMorphy}.
#' @template concavityParam
#' @template pointlessDots
#'
#' @return The 'fit', `h / h + k`, where `h` is the amount of homoplasy ('extra steps') 
#'         and `k` is a constant (the 'concavity constant')
#'
#' @references
#'  - \insertRef{Goloboff1997}{TreeSearch}
#'  
#'  - \insertRef{SmithTern}{TreeSearch}
#'
#' @examples
#'   data(referenceTree)
#'   data(congreveLamsdellMatrices)
#'   dataset <- PrepareDataIW(congreveLamsdellMatrices[[42]])
#'   IWScore(referenceTree, dataset)
#'
#' @author Martin R. Smith
#' @family tree scoring
#' @keywords tree
#' @export
IWScore <- function (tree, dataset, concavity = 10, ...) {
  if (class(dataset) != 'phyDat') {
    stop('Data not of class phyDat; see PhyDat() and PrepareDataIW().')
  }
  if (!('min.length' %in% names(attributes(dataset)))) {
    dataset <- PrepareDataIW(dataset)
  }
  at <- attributes(dataset)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- CharacterLength(tree, dataset)
  minLength <- at$min.length
  homoplasies <- steps - minLength
  
  # This check has been triggered - underlying C memory failure suspected
  # but remains under investigation...
  if (any(homoplasies < 0)) stop("Minimum steps have been miscalculated.\n", 
    "       Please report this bug at:\n", 
    "       https://github.com/ms609/TreeSearch/issues/new\n\n",
    "       See above for full tree: ", dput(tree))
  fit <- homoplasies / (homoplasies + concavity)
  # Return:
  sum(fit * weight)
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
                  parent=parent, child=child, integer(1))
  homoplasies <- steps - minLength
  fit <- homoplasies / (homoplasies + concavity)
  # Return:
  sum(fit * attr(dataset, 'weight'))
}

#' @describeIn IWScore Initialize dataset by adding `morphyObjs` and 
#' `min.length` properties.
#' @export
IWInitMorphy <- function (dataset) {
  attr(dataset, 'morphyObjs') <- 
    lapply(PhyToString(dataset, byTaxon=FALSE, useIndex=FALSE, concatenate=FALSE), 
           SingleCharMorphy)
  
  # Return:
  dataset
}


#' @describeIn TreeSearch Search using profile parsimony
#' @template concavityParam
#' @export
IWTreeSearch <- function (tree, dataset, concavity = 10, 
                          EdgeSwapper = RootedTBR,
                        maxIter = 100, maxHits = 20, forestSize = 1,
                        verbosity = 1, ...) {
  if (class(dataset) != 'phyDat') stop("Unrecognized dataset class; should be phyDat, not ", class(dataset), '.')
  if (!('min.length' %in% names(attributes(dataset)))) dataset <- PrepareDataIW(dataset)
  at <- attributes(dataset)
  
  TreeSearch(tree, dataset, nChar=at$nr, weight=at$weight,
             minLength=at$min.length, concavity = concavity,
             InitializeData = IWInitMorphy,
             CleanUpData = IWDestroyMorphy,
             TreeScorer = IWScoreMorphy,
             EdgeSwapper = EdgeSwapper,
             maxIter = maxIter, maxHits = maxHits, forestSize = forestSize,
             verbosity = verbosity, ...)
}