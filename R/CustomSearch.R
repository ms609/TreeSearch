#' @describeIn TreeSearch Tree search from edge lists
#' @template edgeListParam
#' @template dataForFunction
#' @keywords internal
#' @export
EdgeListSearch <- function (edgeList, dataset,
                            TreeScorer = MorphyLength,
                            EdgeSwapper = RootedTBRSwap,
                            maxIter = 100, maxHits = 20, 
                            bestScore = NULL, stopAtScore = NULL, 
                            stopAtPeak = FALSE, stopAtPlateau = 0L,
                            verbosity = 1L, ...) {
  epsilon <- 1e-07
  
  if (is.null(bestScore)) {
    if (length(edgeList) < 3L) {
      bestScore <- TreeScorer(edgeList[[1]], edgeList[[2]], dataset, ...)
    } else {
      bestScore <- edgeList[[3]]
    }
  }
  if (verbosity > 0L) {
    message("  - Performing tree search.  Initial score: ", bestScore) #nocov
  }
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) {
    if (verbosity > 0L) {  #nocov start
      message("  - Aborting tree search as tree score ", bestScore, 
              " already below target of ", stopAtScore)
    } #nocov end
    edgeList[[3]] <- bestScore
    return(edgeList)
  }
  hits <- 0L
  unimprovedSince <- 0L
  
  for (iter in 1:maxIter) {
    candidateLists <- RearrangeEdges(edgeList[[1]], edgeList[[2]], 
                                     dataset = dataset, 
                                     TreeScorer = TreeScorer,
                                     EdgeSwapper = EdgeSwapper,
                                     hits = hits, iter = iter,
                                     verbosity = verbosity, ...)
    scoreThisIteration <- candidateLists[[3]]
    hits <- candidateLists[[4]]

    if (scoreThisIteration < bestScore + epsilon) {
      if (scoreThisIteration + epsilon < bestScore) unimprovedSince <- -1L
      bestScore <- scoreThisIteration
      edgeList  <- candidateLists
      if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) {
        return(edgeList)
      }
    } else if (stopAtPeak && scoreThisIteration > bestScore + epsilon) {
      if (verbosity > 1L) { #nocov start
        message("    ! Iteration ", iter, 
                " - No TBR rearrangement improves score. ",
                scoreThisIteration, " doesn't beat ", bestScore)
      } #nocov end
      break
    }
    unimprovedSince <- unimprovedSince + 1L
    if (stopAtPlateau > 0L) {
      if (verbosity > 2L && unimprovedSince > 0L) {
        message(" Last improvement ", unimprovedSince, " iterations ago.")
      }
      if (unimprovedSince >= stopAtPlateau) {
        if (verbosity > 1L) message("  - Terminating search, as score has ",
                                    "not improved over past ",
                                    unimprovedSince, " searches.")
        break
      }
    }
    
    if (hits >= maxHits) {
      if (verbosity > 1L) { #nocov start
        message("  - Terminating search; hit best score ", hits, " times.")
      } #nocov end
      break
    }
  }
  if (verbosity > 0L) { #nocov start
    message("  - Final score ", bestScore, " found ", hits, " times after ",
            iter, " rearrangements.", if (verbosity > 1L) '\n' else '')
  } #nocov end
  
  edgeList[3:4] <- c(bestScore, hits)
  
  # Return:
  edgeList
}

#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR})
#' to search for a more parsimonious tree.
#'  
#' For detailed documentation of the 'TreeSearch' package, including full
#' instructions for loading phylogenetic data into R and initiating and
#' configuring tree search, see the
#' [package documentation](https://ms609.github.io/TreeSearch/).
#'  
#' @param tree A fully-resolved starting tree in \code{\link{phylo}} format, 
#' with the desired outgroup.
#' Edge lengths are not supported and will be removed.
#' @template datasetParam
#' @template EdgeSwapperParam
#' @param maxIter Numeric specifying maximum number of iterations to perform
#' before abandoning the search.
#' @param maxHits Numeric specifying maximum times to hit the best pscore
#' before abandoning the search.
#' @template stopAtPeakParam
#' @template stopAtPlateauParam
#'
#' @template InitializeDataParam
#' @template CleanUpDataParam
#' @template treeScorerParam
#'
#' @template verbosityParam
#' @template treeScorerDots
#' 
#' @return
#' `TreeSearch()` returns a tree, with an attribute `pscore` conveying its
#' parsimony score.
#' #' Note that the parsimony score will be inherited from the tree's
#' attributes, which is only valid if it was generated using the same
#' `data` that is passed here.
#' 
#'
#' @references
#' - \insertRef{Smith2019}{TreeSearch}
#'
#' @seealso
#' \itemize{
#' \item \code{\link{Fitch}}, calculates parsimony score;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#' \item \code{\link{Ratchet}}, alternative heuristic, useful to escape local
#'   optima.
#' }
#'
#' @examples
#' data('Lobo', package='TreeTools')
#' njtree <- TreeTools::NJTree(Lobo.phy)
#'
#' \dontrun{
#' TreeSearch(njtree, Lobo.phy, maxIter = 20, EdgeSwapper = NNISwap)
#' TreeSearch(njtree, Lobo.phy, maxIter = 20, EdgeSwapper = RootedSPRSwap)
#' TreeSearch(njtree, Lobo.phy, maxIter = 20, EdgeSwapper = TBRSwap)
#' }
#' @template MRS
#' @family custom search functions
#' @importFrom Rdpack reprompt
#' @importFrom TreeTools RenumberTips
#' @export
TreeSearch <- function (tree, dataset,
                        InitializeData = PhyDat2Morphy,
                        CleanUpData    = UnloadMorphy,
                        TreeScorer     = MorphyLength,
                        EdgeSwapper    = RootedTBRSwap,
                        maxIter = 100L, maxHits = 20L,
                        stopAtPeak = FALSE, stopAtPlateau = 0L,
                        verbosity = 1L, ...) {
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) {
    stop("tree must be bifurcating; try rooting with ape::root")
  }
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- tree$edge
  edgeList <- RenumberEdges(edgeList[, 1], edgeList[, 2])

  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))

  bestScore <- attr(tree, 'score')
  edgeList <- EdgeListSearch(edgeList, initializedData, TreeScorer = TreeScorer, 
                             EdgeSwapper = EdgeSwapper, maxIter = maxIter, 
                             maxHits = maxHits, stopAtPeak = stopAtPeak,
                             stopAtPlateau = stopAtPlateau,
                             verbosity = verbosity, ...)
  
  tree$edge <- cbind(edgeList[[1]], edgeList[[2]])
  attr(tree, 'score') <- edgeList[[3]]
  attr(tree, 'hits') <- edgeList[[4]]
  
  # Return:
  tree 
}
