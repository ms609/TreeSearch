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

#' @describeIn Ratchet Shortcut for Ratchet search using implied weights
#' @template concavityParam
#' @export
IWRatchet <- function (tree, dataset, concavity = 10,
                       swappers = list(TBRSwap, SPRSwap, NNISwap),
                       BootstrapSwapper = if (is.list(swappers))
                         swappers[[length(swappers)]] else swappers,
                       returnAll=FALSE, stopAtScore=NULL,
                       stopAtPeak=FALSE, stopAtPlateau=0L, 
                       ratchIter=100, ratchHits=10, searchIter=2000, searchHits=40,
                       bootstrapIter=searchIter, bootstrapHits=searchHits, verbosity=1L, 
                       suboptimal=1e-08, ...) {
  .Deprecated("MaximizeParsimony")
  dataset <- PrepareDataIW(dataset)
  if (verbosity > 1L) {
    message("* Using implied weighting with concavity constant k = ", concavity)
  }
  
  Ratchet(tree=tree, dataset=dataset, 
          concavity=concavity, minLength=attr(dataset, 'min.length'), 
          InitializeData=IWInitMorphy, CleanUpData=IWDestroyMorphy,
          TreeScorer=IWScoreMorphy, Bootstrapper=IWBootstrap,
          swappers=swappers, BootstrapSwapper=BootstrapSwapper,
          returnAll=returnAll, suboptimal=suboptimal, stopAtScore=stopAtScore,
          ratchIter=ratchIter, ratchHits=ratchHits,
          stopAtPeak=stopAtPeak, stopAtPlateau=stopAtPlateau, 
          searchIter=searchIter, searchHits=searchHits,
          bootstrapIter=searchIter, bootstrapHits=bootstrapHits, 
          verbosity=verbosity, ...)
}

#' @rdname Ratchet 
#' @return `IWMultiRatchet` returns a list of optimal trees produced by `nSearch` 
#'                     Ratchet searches, using implied weighting.
#' @export
IWMultiRatchet <- function (tree, dataset, ratchHits=10, concavity=4,
                            searchIter=500, searchHits=20, verbosity=0L, 
                            swappers=list(RootedNNISwap), nSearch=10, 
                            suboptimal=suboptimal,
                            stopAtScore=NULL, ...) {
  trees <- lapply(seq_len(nSearch), function (i) {
    if (verbosity > 1L) message("\nRatchet search ", i, '/', nSearch, ':')
    IWRatchet(tree, dataset, ratchIter = 1L, ratchHits = 0L, 
              concavity = concavity, 
              searchIter = searchIter, searchHits = searchHits, 
              verbosity = verbosity, swappers = swappers,
              stopAtScore = stopAtScore, ...)
  })
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- UniqueExceptHits(trees[scores == min(scores)])
  message("Found ", length(trees), ' unique trees from ', nSearch, ' searches.')
  
  # Return:
  structure(trees, class = 'multiPhylo')
}