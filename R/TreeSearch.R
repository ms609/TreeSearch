#' @describeIn TreeSearch Tree Search from Edge lists
#' @template edgeListParam
#' @template dataForFunction
#' @author Martin R. Smith
#' @keywords internal
#' @export
EdgeListSearch <- function (edgeList, dataset,
                          TreeScorer = MorphyLength,
                          EdgeSwapper = RootedTBRSwap,
                          maxIter=100, maxHits=20, 
                          bestScore=NULL, stopAtScore=NULL, 
                          stopAtPeak=FALSE, stopAtPlateau=0L,
                          forestSize = 1L, verbosity = 1L, ...) {
  epsilon <- 1e-07
  if (!is.null(forestSize) && length(forestSize)) {
    if (forestSize > 1L) {
      stop("TODO: Forests not supported")
      #### forest <- empty.forest <- vector('list', forestSize)
      #### forest[[1]] <- edgeList
    } else {
      forestSize <- 1L
    }
  }
  if (is.null(bestScore)) {
    if (length(edgeList) < 3L) {
      bestScore <- TreeScorer(edgeList[[1]], edgeList[[2]], dataset, ...)
    } else {
      bestScore <- edgeList[[3]]
    }
  }
  if (verbosity > 0L) {
    message("  - Performing tree search.  Initial score: ", bestScore)
  }
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) {
    if (verbosity > 0L) {
      message("  - Aborting tree search as tree score ", bestScore, 
              " already below target of ", stopAtScore)
    }
    edgeList[[3]] <- bestScore
    return(edgeList)
  }
  returnSingle <- !(forestSize > 1L)
  hits <- 0L
  unimprovedSince <- 0L
  
  for (iter in 1:maxIter) {
    candidateLists <- RearrangeEdges(edgeList[[1]], edgeList[[2]], dataset=dataset, 
                             TreeScorer = TreeScorer, EdgeSwapper=EdgeSwapper, 
                             hits=hits, iter=iter, verbosity=verbosity, ...)
    scoreThisIteration <- candidateLists[[3]]
    hits <- candidateLists[[4]]
    if (forestSize > 1L) {
      stop("TODO re-code this")
      ###if (scoreThisIteration == bestScore) {
      ###  forest[(hits-length(candidateLists)+1L):hits] <- candidateLists ## TODO Check that length still holds
      ###  edgeList  <- sample(forest[1:hits], 1)[[1]]
      ###  bestScore <- scoreThisIteration
      ###  hits      <- hits + 1L
      ###} else if (scoreThisIteration < bestScore) {
      ###  bestScore <- scoreThisIteration
      ###  forest <- empty.forest
      ###  forest[1:hits] <- candidateLists
      ###  edgeList <- sample(candidateLists , 1)[[1]]
      ###  attr(edgeList, 'score') <- scoreThisIteration
      ###}
    } else {
      if (scoreThisIteration < bestScore + epsilon) {
        if (scoreThisIteration + epsilon < bestScore) unimprovedSince <- -1L
        bestScore <- scoreThisIteration
        edgeList  <- candidateLists
        if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(edgeList)
      } else if (stopAtPeak && scoreThisIteration > bestScore + epsilon) {
        if (verbosity > 1L) {
          message("    ! Iteration ", iter, 
                  " - No TBR rearrangement improves score. ",
                  scoreThisIteration, " doesn't beat ", bestScore)
        }
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
    }
    if (hits >= maxHits) {
      if (verbosity > 1L) message("  - Terminating search; hit best score ",
                                  hits, " times.")
      break
    }
  }
  if (verbosity > 0L) {
    message("  - Final score ", bestScore, " found ", hits, " times after ",
            iter, " rearrangements.", if (verbosity > 1L) '\n' else '')
  }
  
  if (forestSize > 1L) {
    if (hits < forestSize) forest <- forest[-((hits+1):forestSize)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- bestScore
    
    # Return:
    unique(forest)
  } else {
    edgeList[3:4] <- c(bestScore, hits)
    
    # Return:
    edgeList
  }
}

#' Search for most parsimonious trees
#'
#' Run standard search algorithms (\acronym{NNI}, \acronym{SPR} or \acronym{TBR}) 
#' to search for a more parsimonious tree.
#'  
#' For detailed documentation of the TreeSearch package, including full 
#' instructions for loading phylogenetic data into R and initiating and 
#' configuring tree search, see the 
#' [package documentation](https://ms609.github.io/TreeSearch).
#'  
#' @param tree a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
#' @template datasetParam
#' @template EdgeSwapperParam
#' @param maxIter the maximum number of iterations to perform before abandoning the search.
#' @param maxHits the maximum times to hit the best pscore before abandoning the search.
#' @template stopAtPeakParam
#' @template stopAtPlateauParam
#' @param forestSize the maximum number of trees to return - useful in concert with \code{\link{consensus}}.
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
#' - \insertRef{SmithTern}{TreeSearch}
#'
#' @seealso
#' \itemize{
#' \item \code{\link{Fitch}}, calculates parsimony score;
#' \item \code{\link{RootedNNI}}, conducts tree rearrangements;
#### \item \code{\link{Sectorial}}, alternative heuristic, useful for larger trees;
#' \item \code{\link{Ratchet}}, alternative heuristic, useful to escape local optima.
#' }
#'
#' @examples
#' data('Lobo', package='TreeTools')
#' njtree <- TreeTools::NJTree(Lobo.phy)
#'
#' \dontrun{
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=NNISwap)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=RootedSPRSwap)
#' TreeSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=TBRSwap)}
#' @template MRS
#' @importFrom Rdpack reprompt
#' @importFrom TreeTools RenumberTips
#' @export
TreeSearch <- function (tree, dataset,
                        InitializeData = PhyDat2Morphy,
                        CleanUpData    = UnloadMorphy,
                        TreeScorer     = MorphyLength,
                        EdgeSwapper    = RootedTBRSwap,
                        maxIter = 100L, maxHits = 20L, forestSize = 1L,
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
                             maxHits = maxHits, forestSize = forestSize, 
                             stopAtPeak = stopAtPeak, stopAtPlateau = stopAtPlateau,
                             verbosity = verbosity, ...)
  
  tree$edge <- cbind(edgeList[[1]], edgeList[[2]])
  attr(tree, 'score') <- edgeList[[3]]
  attr(tree, 'hits') <- edgeList[[4]]
  # Return:
  tree 
}

#' Find most parsimonious trees
#' 
#' Work in progress...
#' 
#' @param ratchHurry Numeric specifying how many times fewer TBR iterations and
#' hits to employ in ratchet searches.
#' @param concavity Numeric specifying concavity constant for implied step 
#' weighting; set as `Inf` for equal step weights (which is a bad idea; see 
#' Smith 2019).
#' 
#'  
#' @examples
#' library(TreeTools)
#' load_all()
#' data('Lobo', package='TreeTools')
#' njtree <- TreeTools::NJTree(Lobo.phy)
#'
#' dataset = Lobo.phy
#' 
#' dataset <- ReadAsPhyDat('c:/research/r/hyoliths/mbank_X24932_6-19-2018_744.nex')
#' tree = NJTree(dataset)
#' 
#' verbosity = 5L
#' tbrIter = 10
#' concavity = 10L
#' session = NULL
#' ratchIter <- tbrIter <- 1L
#' maxHits = 10
#' finalIter = 1
#' 
#' profvis::profvis(MaximizeParsimony(dataset[1:14], concavity = Inf, maxHits = 10, 
#' ratchIter = 3L, tbrIter = 6L, finalIter = 3L))
#' 
#' MaximizeParsimony(dataset, verbosity = 3, concavity = 10, maxHits = 10)
#' 
#' @importFrom TreeTools NJTree
#' @references \insertRef{Smith2019}{TreeTools}
#' @template MRS
#' @export
MaximizeParsimony <- function (dataset, tree = NJTree(dataset),
                               ratchIter = 6L, tbrIter = 6L,
                               maxHits = 100L,
                               finalIter = 3L,
                               concavity = Inf,
                               verbosity = 1L, session = NULL) {
  # Definitions
  .Message <- if (is.null(session)) function (level, ...) {
    if (level < verbosity) {
      message(rep(' ', level), '- ', ...)
    }
  } else function (level, ...) {
    if (level < verbosity) {
      setProgress(message = paste0(...))
    }
  }
  
  .NewOperation <- function(...) if (!is.null(session)) {
    setProgress(0, message = paste0(...))
  }
  .Progress <- function(x, ...) if (!is.null(session)) {
    setProgress(x, message = paste0(...))
  }
  
  .TBRSearch <- function (edge, nTip, morphyObj, tbrIter, maxHits) {
    if (is.null(session)) {
      .DoTBRSearch(edge, nTip, morphyObj, tbrIter, maxHits)
    } else {
      withProgress(message = 'TBR search',
        .DoTBRSearch(edge, nTip, morphyObj, tbrIter, maxHits)
      )
    }
  }
  
  .DoTBRSearch <- function (edge, nTip, morphyObj, tbrIter, maxHits) {
    
    iter <- 0L
    nHits <- 1L
    hold <- array(NA, dim = c(dim(edge), max(maxHits * 1.1, maxHits + 10L)))
    hold[, , 1] <- edge
    bestScore <- preorder_morphy(edge, morphyObj)
    
    while (iter < tbrIter) {
      .Progress(iter / tbrIter, detail = paste0('TBR iteration ', iter + 1))
      iter <- iter + 1L
      optTbr <- sample(3:(nTip * 2 - 2))
      .Message(2L, "TBR iteration ", iter)
      
      for (brk in optTbr) {
        .Message(6L, "Break ", brk)
        moves <- TBRMoves(edge, brk)
        improvedScore <- FALSE
        for (move in moves[sample(seq_along(moves))]) {
          moveScore <- preorder_morphy(move, morphyObj)
          if (moveScore < bestScore + epsilon) {
            edge <- move
            if (moveScore < bestScore) {
              improvedScore <- TRUE
              iter <- 1L
              bestScore <- moveScore
              nHits <- 1L
              hold[, , 1] <- edge
              .Message(1L, "New best score ", bestScore)
              break
            } else {
              .Message(3L, "Best score ", bestScore, " hit again (", nHits, 
                       "/", maxHits, ")")
              nHits <- nHits + 1L
              hold[, , nHits] <- edge
              break
            }
          }
        }
        if (nHits >= maxHits) break
      }
      if (nHits >= maxHits) break
    }
    .Message(0L, "Final score ", bestScore, " found ", nHits, " times after ",
             iter, " rearrangements.")
    
    
    # Return:
    unique(hold[, , seq_len(nHits), drop = FALSE], MARGIN = 3L)
  }
  
  .IWTBRSearch <- function (edge, nTip, morphyObjects, weight, minLength,
                            concavity, tbrIter, maxHits) {
    if (is.null(session)) {
      .DoIWTBRSearch(edge, nTip, morphyObjects, weight, minLength, concavity,
                     tbrIter, maxHits)
    } else {
      withProgress(message = 'TBR search',
                   .DoIWTBRSearch(edge, nTip, morphyObjects, weight, minLength,
                                  concavity, tbrIter, maxHits)
      )
    }
  }
  
  .DoIWTBRSearch <- function (edge, nTip, morphyObjects, weight, minLength,
                              concavity, tbrIter, maxHits) {
    iter <- 0L
    nHits <- 1L
    hold <- array(NA, dim = c(dim(edge), max(maxHits * 1.1, maxHits + 10L)))
    hold[, , 1] <- edge
    bestScore <- .IWScore(edge, morphyObjects, weight, minLength, concavity)
    
    while (iter < tbrIter) {
      .Progress(iter / tbrIter, detail = paste0('TBR iteration ', iter + 1))
      iter <- iter + 1L
      optTbr <- sample(3:(nTip * 2 - 2))
      .Message(2L, "TBR iteration ", iter)
      
      for (brk in optTbr) {
        .Message(6L, "Break ", brk)
        moves <- TBRMoves(edge, brk)
        improvedScore <- FALSE
        for (move in moves[sample(seq_along(moves))]) {
          moveScore <- .IWScore(move, morphyObjects, weight, minLength, concavity)
          if (moveScore < bestScore + epsilon) {
            edge <- move
            if (moveScore < bestScore) {
              improvedScore <- TRUE
              iter <- 0L
              bestScore <- moveScore
              nHits <- 1L
              hold[, , 1] <- edge
              .Message(1L, "New best score ", signif(bestScore, 5))
              break
            } else {
              .Message(3L, "Best score ", signif(bestScore, 5),
                       " hit again (", nHits, "/", maxHits, ")")
              nHits <- nHits + 1L
              hold[, , nHits] <- edge
              break
            }
          }
        }
        if (nHits >= maxHits) break
      }
      if (nHits >= maxHits) break
    }
    .Message(0L, "Final score ", signif(bestScore, 5), " found ", 
             nHits, " times after ", iter, " rearrangements.")
    
    
    # Return:
    unique(hold[, , seq_len(nHits), drop = FALSE], MARGIN = 3L)
  }
  
  .IWScore <- function (edge, morphyObjs, weight, minLength, concavity) {
    steps <- preorder_morphy_by_char(edge, morphyObjs)
    homoplasies <- steps - minLength
    fit <- homoplasies / (homoplasies + concavity)
    sum(fit * weight)
  }
  
  epsilon <- sqrt(.Machine$double.eps)
  iw <- is.finite(concavity)
  
  # Initialize tree
  if (inherits(tree, 'multiPhylo')) {
    .Message(1L, "Starting search from `tree[[1]]`.")
    tree <- tree[[1]]
  }
  if (dim(tree$edge)[1] != 2 * tree$Nnode) {
    stop("`tree` must be bifurcating; try rooting with RootTree(tree, root = 1)")
  }
  
  tree <- Preorder(RenumberTips(tree, names(dataset)))
  nTip <- NTip(tree)
  edge <- tree$edge
  
  # Initialize data
  if (iw) {
    at <- attributes(dataset)
    characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                              byTaxon = FALSE, concatenate = FALSE)
    startWeights <- at$weight
    morphyObjects <- lapply(characters, SingleCharMorphy)
    on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)))
    
    nLevel <- length(at$level)
    nChar <- at$nr
    nTip <- length(dataset)
    cont <- at$contrast
    inappLevel <- at$levels == '-'
    
    if (any(inappLevel)) {
      # TODO this is a workaround until MinimumLength can handle {-, 1}
      cont[cont[, inappLevel] > 0, ] <- 0
      ambiguousToken <- at$allLevels == '?'
      cont[ambiguousToken, ] <- colSums(cont[!ambiguousToken, ]) > 0
    }
    
    # Perhaps replace with previous code:
    # inappLevel <- which(at$levels == "-")
    # cont[, inappLevel] <- 0
    
    powersOf2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
    tmp <- as.integer(cont %*% powersOf2)
    unlisted <- unlist(dataset, use.names = FALSE)
    binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow = FALSE)
    
    minLength <- apply(binaryMatrix, 1, MinimumLength)
  } else {
    morphyObj <- PhyDat2Morphy(dataset)
    on.exit(morphyObj <- UnloadMorphy(morphyObj))
    startWeights <- MorphyWeights(morphyObj)[1, ]
  }
  
  # Prepare search
  .Message(0L, "Parsimony search with ", ratchIter, " ratchet iterations; ", 
           tbrIter, " TBR rounds; ", maxHits, " hits; k = ", concavity, ".")
  
  iter <- 0L
  nHits <- 1L
  if (ratchIter > 0L) {
    .Message(0L, "Performing parsimony ratchet.")
    verbosity <- verbosity - 1L
    while (iter < ratchIter) {
      .Progress(iter / ratchIter, "Ratchet iteration ", (iter + 1L))
      iter <- iter + 1L
      .Message(1L, "Ratchet iteration ", iter)
      eachChar <- seq_along(startWeights)
      deindexedChars <- rep(eachChar, startWeights)
      resampling <- tabulate(sample(deindexedChars, replace = TRUE),
                             length(startWeights))
      if (iw) {
        sampled <- resampling > 0L
        ratchetObjs <- morphyObjects[sampled]
        ratchetTrees <- .IWTBRSearch(edge, NTip(tree), ratchetObjs,
                                     resampling[sampled], minLength[sampled],
                                     concavity, tbrIter, maxHits / finalIter)
      } else {
        errors <- vapply(eachChar, function (i) 
          mpl_set_charac_weight(i, resampling[i], morphyObj), integer(1))
        if (any(errors)) {
          stop ("Error resampling morphy object: ",
                mpl_translate_error(unique(errors[errors < 0L])))
        }
        if (mpl_apply_tipdata(morphyObj) -> error) {
          stop("Error applying tip data: ", mpl_translate_error(error))
        }
        
        ratchetTrees <- .TBRSearch(edge, NTip(tree), morphyObj, tbrIter,
                                   maxHits / finalIter)
        errors <- vapply(eachChar, function (i) 
          mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
        if (any(errors)) stop ("Error resampling morphy object: ",
                               mpl_translate_error(unique(errors[errors < 0L])))
        if (mpl_apply_tipdata(morphyObj) -> error) {
          stop("Error applying tip data: ", mpl_translate_error(error))
        }
      }
    }
    verbosity <- verbosity + 1L
    edge <- ratchetTrees[, , sample.int(dim(ratchetTrees)[3], 1)]
    bestEdges <- if (iw) {
      .IWTBRSearch(edge, NTip(tree), morphyObjects, startWeights, minLength,
                   concavity, tbrIter, maxHits / finalIter)
    } else {
      .TBRSearch(edge, NTip(tree), morphyObj, tbrIter, maxHits / finalIter)
    }
    edge <- bestEdges[, , sample.int(dim(bestEdges)[3], 1)]
  }
  
  # Branch breaking
  
  .Message(0L, "Final TBR search.")
  bestEdges <- if (iw) {
    .IWTBRSearch(edge, nTip, morphyObjects, startWeights, minLength,
                 concavity, tbrIter * finalIter, maxHits)
  } else {
    .TBRSearch(edge, nTip, morphyObj, tbrIter * finalIter, maxHits)
  }
  
  ret <- structure(lapply(seq_len(dim(bestEdges)[3]), function (i) {
    tr <- tree
    tr$edge <- bestEdges[, , i]
    tr
  }), class = 'multiPhylo')
  
  # Return:
  ret
}
