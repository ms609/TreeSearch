#' Find most parsimonious trees
#' 
#' Search for most parsimonious trees using the parsimony ratchet and 
#' \acronym{TBR} rearrangements, treating inapplicable data as such using the
#' algorithm of Brazeau, Guillerme & Smith (2019).
#'  
#' Tree search will be conducted from a specified or automatically-generated
#' starting tree in order to find a tree with an optimal parsimony score,
#' under implied or equal weights, treating inapplicable characters as such
#' in order to avoid the artefacts of the standard Fitch algorithm
#' (see Maddison 1993; Brazeau et al. 2019).
#' 
#' Tree search commences with `ratchIter` iterations of the parsimony ratchet
#' (Nixon 1999), which bootstraps the input dataset in order to escape local
#' optima.  A final round of tree bisection and reconnection (\acronym{TBR})
#' is conducted to broaden the sampling of trees.
#' 
#' This function can be called using the R command line / terminal, or through
#' the 'shiny' graphical user interface app (type `EasyTrees()` to launch).
#' 
#'  
#' For detailed documentation of the TreeSearch package, including full 
#' instructions for loading phylogenetic data into R and initiating and 
#' configuring tree search, see the 
#' [package documentation](https://ms609.github.io/TreeSearch).
#'  
#' 
#' @template datasetParam
#' @param tree (optional) A bifurcating tree of class \code{\link{phylo}},
#' containing only the tips listed in `dataset`, from which the search
#' should begin.
#' If unspecified, a neighbour-joining tree will be generated from `dataset`.
#' Edge lengths are not supported and will be deleted.
#' @param ratchIter Numeric specifying number of iterations of the 
#' parsimony ratchet (Nixon 1999) to conduct.
#' @param tbrIter Numeric specifying the maximum number of \acronym{TBR}
#' break points to evaluate before concluding each search.
#' The counter is reset to zero each time tree score improves.
#' One 'iteration' comprises breaking a single branch and evaluating all 
#' possible reconnections.
#' @param finalIter Numeric: the final round of tree search will evaluate
#' `finalIter` &times; `tbrIter` \acronym{TBR} break points.
#' @param maxHits Numeric specifying the maximum times that an optimal
#' parsimony score may be hit before concluding a ratchet iteration or final 
#' search concluded.
#' @param concavity Numeric specifying concavity constant for implied step 
#' weighting; set as `Inf` for equal step weights (which is a bad idea; see 
#' Smith 2019).
#' @param tolerance Numeric specifying degree of suboptimality to tolerate
#' before rejecting a tree.  The default, `sqrt(.Machine$double.eps)`, retains
#' trees that may be equally parsimonious but for rounding errors.  
#' Setting to larger values will include trees suboptimal by up to `tolerance`
#' in search results, which may improve the accuracy of the consensus tree
#' (at the expense of resolution) (Smith 2019).
#' @param constraint Either `NULL` or an object of class `phyDat`. Trees that
#' are not perfectly compatible with each character in `constraint` will not
#' be considered during search.
#' See [vignette](https://ms609.github.io/TreeSearch/articles/inapplicable.html)
#' for further examples.
#' @param verbosity Integer specifying level of messaging; higher values give
#' more detailed commentary on search progress. Set to `0` to run silently.
#' @param session 'shiny' session identifier to allow [`setProgress()`] calls
#' to be sent when `MaximizeParsimony()` is called from within a shiny app..
#' 
#'  
#' @examples
#' library('TreeTools')
#' data('Lobo', package = 'TreeTools')
#' dataset <- Lobo.phy
#' 
#' # A very quick run for demonstration purposes
#' MaximizeParsimony(dataset, ratchIter = 0, tbrIter = 1, concavity = 10,
#'                   maxHits = 5, verbosity = 4)
#' 
#' # Be sure to check that the score has converged on a global optimum,
#' # conducting additional iterations and runs as necessary.
#' 
#' \dontrun{ # launches 'shiny' point-and-click interface
#'   EasyTrees() 
#' }
#' 
#' # Tree search with a constraint
#' constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
#' characters <- MatrixToPhyDat(matrix(
#'   c(0, 1, 1, 1, 0, 0,
#'     1, 1, 1, 0, 0, 0), ncol = 2,
#'   dimnames = list(letters[1:6], NULL)))
#' MaximizeParsimony(characters, constraint = constraint, verbosity = 0)
#' 
#' @template MRS
#' 
#' @importFrom phangorn Descendants
#' @importFrom shiny setProgress withProgress
#' @importFrom TreeTools NJTree CharacterInformation
#' @references
#' \insertRef{Brazeau2019}{TreeSearch}
#' 
#' \insertRef{Maddison1993}{TreeSearch}
#' 
#' \insertRef{Nixon1999}{TreeSearch}
#' 
#' \insertRef{Smith2019}{TreeSearch}
#' @export
MaximizeParsimony <- function (dataset, tree = NJTree(dataset),
                               ratchIter = 12L, tbrIter = 6L, finalIter = 3L,
                               maxHits = 20L,
                               concavity = Inf,
                               tolerance = sqrt(.Machine$double.eps),
                               constraint = NULL,
                               verbosity = 2L, session = NULL) {
  # Define functions
  .Message <- if (is.null(session)) function (level, ...) {
    if (level < verbosity) {
      message(rep(' ', level), '- ', ...)
    }
  } else function (level, ...) { # nocov start
    if (level < verbosity) {
      setProgress(message = paste0(...))
    }
  } # nocov end
  
  .NewOperation <- function(...) if (!is.null(session)) { # nocov start
    setProgress(0, message = paste0(...))
  }
  .Progress <- function(x, ...) if (!is.null(session)) {
    setProgress(x, message = paste0(...))
  } # nocov end
  
  .TBRSearch <- function (edge, nTip, morphyObj, tbrIter, maxHits) {
    if (is.null(session)) {
      .DoTBRSearch(edge, nTip, morphyObj, tbrIter, maxHits)
    } else { # nocov start
      withProgress(message = 'TBR search',
                   .DoTBRSearch(edge, nTip, morphyObj, tbrIter, maxHits)
      )
    } # nocov end
  }
  
  .DoTBRSearch <- function (edge, nTip, morphyObj, tbrIter, maxHits) {
    
    iter <- 0L
    nHits <- 1L
    hold <- array(NA, dim = c(dim(edge), max(maxHits * 1.1, maxHits + 10L)))
    hold[, , 1] <- edge
    bestScore <- preorder_morphy(edge, morphyObj)
    
    while (iter < tbrIter) {
      .Progress(iter / tbrIter, detail = paste0('TBR iteration (depth ', iter + 1, ')'))
      iter <- iter + 1L
      optTbr <- sample(3:(nTip * 2 - 2))
      .Message(3L, "New TBR iteration (depth ", iter, ')')
      
      for (brk in optTbr) {
        .Message(6L, "Break ", brk)
        moves <- TBRMoves(edge, brk)
        improvedScore <- FALSE
        for (move in moves[sample(seq_along(moves))]) {
          if (.Forbidden(move)) {
            .Message(10L, "Skipping prohibited topology")
            next
          }
          moveScore <- preorder_morphy(move, morphyObj)
          if (moveScore < bestScore + epsilon) {
            edge <- move
            if (moveScore < bestScore) {
              improvedScore <- TRUE
              iter <- 1L
              bestScore <- moveScore
              nHits <- 1L
              hold[, , 1] <- edge
              .Message(4L, "New best score ", bestScore)
              break
            } else {
              .Message(5L, "Best score ", bestScore, " hit again (", nHits, 
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
    .Message(2L, iter + 1, " TBR rearrangements found score ",
             signif(bestScore), " ", nHits, " times.")
    
    
    # Return:
    unique(hold[, , seq_len(nHits), drop = FALSE], MARGIN = 3L)
  }
  
  .IWTBRSearch <- function (edge, nTip, morphyObjects, weight, minLength,
                            charSeq, concavity, tbrIter, maxHits) {
    if (is.null(session)) {
      .DoIWTBRSearch(edge, nTip, morphyObjects, weight, minLength, charSeq,
                     concavity, tbrIter, maxHits)
    } else { # nocov start
      withProgress(message = 'TBR search',
                   .DoIWTBRSearch(edge, nTip, morphyObjects, weight, minLength,
                                  charSeq, concavity, tbrIter, maxHits)
      )
    } # nocov end
  }
  
  .DoIWTBRSearch <- function (edge, nTip, morphyObjects, weight, minLength,
                              charSeq, concavity, tbrIter, maxHits) {
    iter <- 0L
    nHits <- 1L
    hold <- array(NA, dim = c(dim(edge), max(maxHits * 1.1, maxHits + 10L)))
    maxHits <- ceiling(maxHits)
    hold[, , 1] <- edge
    bestScore <- .IWScore(edge, morphyObjects, weight, minLength, charSeq,
                          concavity)
    bestPlusEps <- bestScore + epsilon
    
    while (iter < tbrIter) {
      .Progress(iter / tbrIter, detail = paste0('TBR iteration (depth ', iter + 1, ')'))
      iter <- iter + 1L
      optTbr <- sample(3:(nTip * 2 - 2))
      .Message(3L, "New TBR iteration (depth ", iter, ')')
      
      for (brk in optTbr) {
        .Message(6L, "Break ", brk)
        moves <- TBRMoves(edge, brk)
        improvedScore <- FALSE
        for (move in moves[sample(seq_along(moves))]) {
          if (.Forbidden(move)) {
            .Message(10L, "Skipping prohibited topology")
            next
          }
          moveScore <- .IWScore(move, morphyObjects, weight, minLength, charSeq,
                                concavity, bestPlusEps)
          if (moveScore < bestPlusEps) {
            edge <- move
            if (moveScore < bestScore) {
              improvedScore <- TRUE
              iter <- 0L
              bestScore <- moveScore
              bestPlusEps <- bestScore + epsilon
              nHits <- 1L
              hold[, , 1] <- edge
              .Message(4L, "New best score ", signif(bestScore, 5))
              break
            } else {
              .Message(5L, "Best score ", signif(bestScore, 5),
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
    .Message(2L, iter + 1, " TBR rearrangements found score ", 
             signif(bestScore), " ", nHits, " times.")
    
    # Return:
    unique(hold[, , seq_len(nHits), drop = FALSE], MARGIN = 3L)
  }
  
  .ProfileTBRSearch <- function (edge, nTip, morphyObjects, weight, charSeq,
                                 profiles, tbrIter, maxHits) {
    if (is.null(session)) {
      .DoProfileTBRSearch(edge, nTip, morphyObjects, weight, charSeq, profiles,
                          tbrIter, maxHits)
    } else { # nocov start
      withProgress(message = 'TBR search',
                   .DoProfileTBRSearch(edge, nTip, morphyObjects, weight, 
                                       charSeq, profiles, tbrIter, maxHits)
      )
    } # nocov end
  }
  
  .DoProfileTBRSearch <- function (edge, nTip, morphyObjects, weight, charSeq,
                                   profiles, tbrIter, maxHits) {
    iter <- 0L
    nHits <- 1L
    hold <- array(NA, dim = c(dim(edge), max(maxHits * 1.1, maxHits + 10L)))
    maxHits <- ceiling(maxHits)
    hold[, , 1] <- edge
    bestScore <- .ProfileScore(edge, morphyObjects, weight, charSeq, profiles)
    bestPlusEps <- bestScore + epsilon
    
    while (iter < tbrIter) {
      .Progress(iter / tbrIter, detail = paste0('TBR iteration (depth ',
                                                iter + 1, ')'))
      iter <- iter + 1L
      optTbr <- sample(3:(nTip * 2 - 2))
      .Message(3L, "New TBR iteration (depth ", iter, ')')
      
      for (brk in optTbr) {
        .Message(6L, "Break ", brk)
        moves <- TBRMoves(edge, brk)
        improvedScore <- FALSE
        for (move in moves[sample(seq_along(moves))]) {
          if (.Forbidden(move)) {
            .Message(10L, "Skipping prohibited topology")
            next
          }
          moveScore <- .ProfileScore(move, morphyObjects, weight, charSeq,
                                     profiles, bestPlusEps)
          if (moveScore < bestPlusEps) {
            edge <- move
            if (moveScore < bestScore) {
              improvedScore <- TRUE
              iter <- 0L
              bestScore <- moveScore
              bestPlusEps <- bestScore + epsilon
              nHits <- 1L
              hold[, , 1] <- edge
              .Message(4L, "New best score ", signif(bestScore, 5))
              break
            } else {
              .Message(5L, "Best score ", signif(bestScore, 5),
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
    .Message(2L, iter + 1, " TBR rearrangements found score ", 
             signif(bestScore), " ", nHits, " times.")
    
    # Return:
    unique(hold[, , seq_len(nHits), drop = FALSE], MARGIN = 3L)
  }
  
  .IWScore <- function (edge, morphyObjs, weight, minLength, charSeq,
                        concavity, target = Inf) {
    morphy_iw(edge, morphyObjs, weight, minLength, charSeq,
              concavity, target + epsilon)
  } 
  
  .ProfileScore <- function (edge, morphyObjs, weight, charSeq, profiles, 
                             target = Inf) {
    morphy_profile(edge, morphyObjs, weight, charSeq, profiles,
                   target + epsilon)
  }
  
  # Initialize tree
  if (inherits(tree, 'multiPhylo')) {
    .Message(1L, "Starting search from `tree[[1]]`.")
    tree <- tree[[1]]
  }
  if (dim(tree$edge)[1] != 2 * tree$Nnode) {
    stop("`tree` must be bifurcating; try rooting with RootTree(tree, 1)")
  }
  
  tree <- Preorder(RenumberTips(tree, names(dataset)))
  nTip <- NTip(tree)
  edge <- tree$edge
  
  if (edge[1, 2] > nTip) {
    outgroup <- Descendants(tree, edge[1, 2], type = 'tips')[[1]]
    if (length(outgroup) > nTip / 2L) {
      outgroup <- seq_len(nTip)[-outgroup]
    }
    tree <- RootTree(tree, 1)
    edge <- tree$edge
  } else {
    outgroup <- NA
  }
  
  # Define constants
  epsilon <- tolerance #sqrt(.Machine$double.eps)
  profile <- tolower(concavity) == "profile"
  iw <- is.finite(concavity)
  constrained <- !is.null(constraint)
  
  # Initialize constraints
  if (constrained) {
    morphyConstr <- PhyDat2Morphy(constraint)
    on.exit(morphyConstr <- UnloadMorphy(morphyConstr), add = TRUE)
    # Calculate constraint minimum score
    constraintLength <- sum(MinimumLength(constraint))
    
    .Forbidden <- function (edges) {
      preorder_morphy(edges, morphyConstr) != constraintLength
    }
    
    # Check that starting tree is consistent with constraints 
    if (.Forbidden(edge)) {
      .Message(1L, "Looking for a tree that is consistent with `constraint`...")
      .oldMsg <- .Message
      .Message <- function (...) {}
      edge <- .DoTBRSearch(edge, nTip, morphyConstr, 10, 10)[, , 1]
      .Message <- .oldMsg
      if (.Forbidden(edge)) {
        stop("Specify a starting tree that is consistent with `constraint`.")
      }
    }
    
  } else {
    .Forbidden <- function (edges) FALSE
  }
  
  # Initialize data
  if (profile) {
    dataset <- PrepareDataProfile(dataset)
    originalLevels <- attr(dataset, 'levels')
    if ('-' %in% originalLevels) {
      #TODO Fixing this will require updating the counts table cleverly
      # Or we could use approximate info amounts, e.g. by treating '-' as 
      # an extra token
      message("Inapplicable tokens '-' treated as ambiguous '?' for profile parsimony")
      cont <- attr(dataset, 'contrast')
      cont[cont[, '-'] != 0, ] <- 1
      attr(dataset, 'contrast') <- cont[, colnames(cont) != '-']
      attr(dataset, 'levels') <- originalLevels[originalLevels != '-']
    }
    profiles <- attr(dataset, 'info.amounts')
  }
  if (iw || profile) {
    at <- attributes(dataset)
    characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                              byTaxon = FALSE, concatenate = FALSE)
    startWeights <- at$weight
    morphyObjects <- lapply(characters, SingleCharMorphy)
    on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)),
            add = TRUE)
    
    nLevel <- length(at$level)
    nChar <- at$nr
    nTip <- length(dataset)
    cont <- at$contrast
    if (is.null(colnames(cont))) colnames(cont) <- as.character(at$levels)
    simpleCont <- ifelse(rowSums(cont) == 1,
                         apply(cont != 0, 1, function (x) colnames(cont)[x][1]),
                         '?')
  }
  if (iw) {
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
  }
  
  if (iw || profile) {
    powersOf2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
    tmp <- as.integer(cont %*% powersOf2)
    unlisted <- unlist(dataset, use.names = FALSE)
    binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow = FALSE)
    minLength <- apply(binaryMatrix, 1, MinimumLength)
    
    tokenMatrix <- matrix(simpleCont[unlisted], nChar, nTip, byrow = FALSE)
    charInfo <- apply(tokenMatrix, 1, CharacterInformation)
    needsInapp <- rowSums(tokenMatrix == '-') > 2
    inappSlowdown <- 3L # A guess
    # Crude estimate of score added per unit processing time
    rawPriority <- charInfo / ifelse(needsInapp, inappSlowdown, 1)
    priority <- startWeights * rawPriority
    informative <- needsInapp | charInfo > 0
    # Will work from end of sequence to start.
    charSeq <- seq_along(charInfo)[informative][order(priority[informative])] - 1L
  } else {
    morphyObj <- PhyDat2Morphy(dataset)
    on.exit(morphyObj <- UnloadMorphy(morphyObj), add = TRUE)
    startWeights <- unlist(MorphyWeights(morphyObj)[1, ]) # exact == approx
  }
  
  # Initialize variables and prepare search
  
  iter <- 0L
  nHits <- 1L
  bestScore <- if (profile) {
    .ProfileScore(edge, morphyObjects, startWeights, charSeq, profiles)
  } else if (iw) {
    .IWScore(edge, morphyObjects, startWeights, minLength, charSeq, concavity)
  } else {
    preorder_morphy(edge, morphyObj)
  }
  bestPlusEps <- bestScore + epsilon
  
  .Message(0L, "Parsimony search with ", ratchIter, " ratchet iterations; ", 
           "TBR depth ", tbrIter, "; ", maxHits, " hits; k = ", concavity, ".",
           "\n  ", Sys.time(),
           "\n  Initial score: ", signif(bestScore, 5),
           "\n")
  if (ratchIter > 0L) {
    .Message(0L, "Performing parsimony ratchet.")
    while (iter < ratchIter) {
      .Progress(iter / ratchIter, "Ratchet iteration ", (iter + 1L))
      iter <- iter + 1L
      .Message(1L, "Ratchet iteration ", iter, 
               ": Search from bootstrapped dataset.")
      verbosity <- verbosity - 1L
      eachChar <- seq_along(startWeights)
      deindexedChars <- rep(eachChar, startWeights)
      resampling <- tabulate(sample(deindexedChars, replace = TRUE),
                             length(startWeights))
      if (profile || iw) {
        priority <- resampling * rawPriority
        sampled <- informative & resampling > 0
        ratchSeq <- seq_along(charInfo)[sampled][order(priority[sampled])] - 1L
      }
      if (profile) {
        ratchetTrees <- .ProfileTBRSearch(edge, NTip(tree), morphyObjects,
                                          resampling, ratchSeq, profiles,
                                          tbrIter, maxHits / finalIter)
      } else if (iw) {
        ratchetTrees <- .IWTBRSearch(edge, NTip(tree), morphyObjects,
                                     resampling, minLength, ratchSeq,
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
      verbosity <- verbosity + 1L
      ratchetStart <- ratchetTrees[, , sample.int(dim(ratchetTrees)[3], 1)]
      .Message(1L, "Ratchet iteration ", iter, ": Search using original data")
      ratchetImproved <- if (profile) {
        .ProfileTBRSearch(ratchetStart, NTip(tree), morphyObjects, startWeights,
                          charSeq, profiles, tbrIter, maxHits / finalIter)
      } else if (iw) {
        .IWTBRSearch(ratchetStart, NTip(tree), morphyObjects, startWeights,
                     minLength, charSeq,
                     concavity, tbrIter, maxHits / finalIter)
      } else {
        .TBRSearch(ratchetStart, NTip(tree), morphyObj, tbrIter, maxHits / finalIter)
      }
      ratchetScore <- if (profile) {
        .ProfileScore(ratchetImproved[, , 1], morphyObjects, startWeights,
                      charSeq, profiles, bestPlusEps)
        
      } else if (iw) {
        .IWScore(ratchetImproved[, , 1], morphyObjects, startWeights, minLength,
                 charSeq, concavity, bestPlusEps)
      } else {
        preorder_morphy(ratchetImproved[, , 1], morphyObj)
      }
      if (ratchetScore < bestPlusEps) {
        if (ratchetScore + epsilon < bestScore) {
          .Message(1L, "*Ratchet iteration ", iter, " found new best score: ",
                   signif(ratchetScore, 5), "*")
          bestScore <- ratchetScore
          bestPlusEps <- bestScore + epsilon
          edge <- ratchetImproved[, , sample.int(dim(ratchetImproved)[3], 1)]
        } else {
          .Message(1L, "Ratchet iteration ", iter, " hit best score ",
                   signif(bestScore, 5), ' again')
          #TODO add ratchet best trees to edge library
          edge <- ratchetImproved[, , sample.int(dim(ratchetImproved)[3], 1)]
        }
      } else {
        .Message(1L, "Ratchet iteration ", iter, " did not hit best score ",
                 signif(bestScore, 5), ".")
      }
    }
  }
  
  # Branch breaking
  
  .Message(0L, "Final TBR search.")
  bestEdges <- if (profile) {
    .ProfileTBRSearch(edge, nTip, morphyObjects, startWeights, charSeq, 
                      profiles, tbrIter * finalIter, maxHits)
  } else if (iw) {
    .IWTBRSearch(edge, nTip, morphyObjects, startWeights, minLength, charSeq,
                 concavity, tbrIter * finalIter, maxHits)
  } else {
    .TBRSearch(edge, nTip, morphyObj, tbrIter * finalIter, maxHits)
  }
  
  finalScore <- if (profile) {
    .ProfileScore(bestEdges[, , 1], morphyObjects, startWeights, charSeq,
                  profiles)
  } else if (iw) {
    .IWScore(bestEdges[, , 1], morphyObjects, startWeights, minLength, charSeq,
             concavity)
  } else {
    preorder_morphy(bestEdges[, , 1], morphyObj)
  }
  .Message(0L, "Final score: ", finalScore, "\n\n")
  
  
  ret <- structure(lapply(seq_len(dim(bestEdges)[3]), function (i) {
    tr <- tree
    tr$edge <- bestEdges[, , i]
    if (is.na(outgroup)) {
      tr
    } else {
      RootTree(tr, outgroup)
    }
  }), class = 'multiPhylo')
  
  # Return:
  ret
}


#' @rdname MaximizeParsimony
Resample <- function (dataset, tree = NJTree(dataset), method = 'jack',
                      ratchIter = 1L, tbrIter = 6L, finalIter = 3L,
                      maxHits = 10L, concavity = Inf,
                      tolerance = sqrt(.Machine$double.eps),
                      constraint = NULL,
                      verbosity = 2L, session = NULL) {
  
}

#' Launch tree search graphical user interface
#' 
#' @rdname MaximizeParsimony
#' @importFrom shiny runApp
#' @importFrom shinyjs useShinyjs
## @importFrom rgl plot3d
#' @importFrom TreeDist ClusteringInfoDistance
#' @importFrom protoclust protoclust
#' @importFrom cluster pam silhouette
#' @export
EasyTrees <- function () # nocov start
  shiny::runApp(system.file('Parsimony', package = 'TreeSearch'))

#' @rdname MaximizeParsimony
#' @export
EasyTreesy <- EasyTrees
# nocov end
