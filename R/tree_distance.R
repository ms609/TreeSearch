#' Information-based generalized Robinson-Fould distance between two trees
#'
#' @param tree1,tree2 Trees of class `phylo`, with tips labelled identically,
#' or lists of such trees to undergo pairwise comparison.
#' 
#' @param reportMatching Logical specifying whether to return the clade
#' matchings as an attribute of the score.
#'
#' @return Returns a numeric that sums the mutual information content of the
#' optimal matching of bipartitions between two trees, following Smith (submitted).
#' 
#' @references {
#'  \insertRef{SmithDist}{TreeSearch}
#' }
#' 
#' @author Martin R. Smith
#' 
#' @importFrom clue solve_LSAP
#' @export
MutualArborealInfo <- function (tree1, tree2, reportMatching = FALSE) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      MutualArborealInfoSplits(Tree2Splits(tree1), Tree2Splits(tree2), reportMatching)
    } else {
      splits1 <- Tree2Splits(tree1)
      vapply(tree2, 
             function (tr2) MutualArborealInfoSplits(splits1, Tree2Splits(tr2)),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- Tree2Splits(tree2)
      vapply(tree1, 
             function (tr2) MutualArborealInfoSplits(splits1, Tree2Splits(tr2)),
             double(1))
    } else {
      splits1 <- lapply(tree1, Tree2Splits)
      splits2 <- lapply(tree2, Tree2Splits)
      matrix(mapply(MutualArborealInfoSplits, rep(splits1, each=length(splits2)), splits2),
             length(splits2), length(splits1), dimnames = list(names(tree2), names(tree1)))
    }
  }
}

#' Variation of phylogenetic information between two trees
#' 
#' See Meila 2007... #TODO document
#' 
#' NOTE: We're using phylogenetic information rather than partitioning
#' information, which imposes arboreal matching.
#' 
#' @references {
#'   \insertRef{Meila2007}{TreeSearch}
#' }
#'
#' @author Martin R. Smith
#' @inheritParams MutualArborealInfo
#' @export
VariationOfArborealInfo <- function (tree1, tree2, reportMatching = FALSE) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      VariationOfSplitArborealInfo(Tree2Splits(tree1), Tree2Splits(tree2), reportMatching)
    } else {
      splits1 <- Tree2Splits(tree1)
      vapply(tree2, 
             function (tr2) VariationOfSplitArborealInfo(splits1, Tree2Splits(tr2)),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- Tree2Splits(tree2)
      vapply(tree1, 
             function (tr2) VariationOfSplitArborealInfo(splits1, Tree2Splits(tr2)),
             double(1))
    } else {
      splits1 <- lapply(tree1, Tree2Splits)
      splits2 <- lapply(tree2, Tree2Splits)
      matrix(mapply(VariationOfSplitArborealInfo, rep(splits1, each=length(splits2)), splits2),
             length(splits2), length(splits1), dimnames = list(names(tree2), names(tree1)))
    }
  }
}

#' Variation of information between two trees
#' 
#' See Meila 2007... #TODO document
#' 
#' NOTE: We're using partition information rather than phylgoenetic information,
#' which does NOT impose arboreal matching.
#'
#' 
#' 
#' @references {
#'   \insertRef{Meila2007}{TreeSearch}
#' }
#'
#' @author Martin R. Smith
#' @inheritParams MutualArborealInfo
#' @export
VariationOfPartitionInfo <- function (tree1, tree2,
                                      reportMatching = FALSE,
                                      bestMatchOnly = TRUE) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      VariationOfSplitPartitionInfo(Tree2Splits(tree1), Tree2Splits(tree2),
                                    reportMatching, bestMatchOnly)
    } else {
      splits1 <- Tree2Splits(tree1)
      vapply(tree2, 
             function (tr2) VariationOfSplitPartitionInfo(splits1, Tree2Splits(tr2),
                                                          reportMatching,
                                                          bestMatchOnly),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- Tree2Splits(tree2)
      vapply(tree1, 
             function (tr2) VariationOfSplitPartitionInfo(splits1, Tree2Splits(tr2),
                                                          reportMatching, bestMatchOnly),
             double(1))
    } else {
      splits1 <- lapply(tree1, Tree2Splits)
      splits2 <- lapply(tree2, Tree2Splits)
      matrix(mapply(VariationOfSplitPartitionInfo,
                    rep(splits1, each=length(splits2)),
                    splits2,
                    reportMatching=reportMatching,
                    bestMatchOnly= bestMatchOnly),
             length(splits2), length(splits1), dimnames = list(names(tree2), names(tree1)))
    }
  }
}

#' Tree distance based on joint information content of splits
#' 
#' #TODO Needs documenting and describing fully.
#' 
#' @inheritParams MutualArborealInfo
#' @author Martin R. Smith
#' @export
MutualPartitionInfo <- function (tree1, tree2, reportMatching = FALSE) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      MutualPartitionInfoSplits(Tree2Splits(tree1), Tree2Splits(tree2), reportMatching)
    } else {
      splits1 <- Tree2Splits(tree1)
      vapply(tree2, 
             function (tr2) MutualPartitionInfoSplits(splits1, Tree2Splits(tr2), reportMatching),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- Tree2Splits(tree2)
      vapply(tree1, 
             function (tr2) MutualPartitionInfoSplits(splits1, Tree2Splits(tr2), reportMatching),
             double(1))
    } else {
      splits1 <- lapply(tree1, Tree2Splits)
      splits2 <- lapply(tree2, Tree2Splits)
      matrix(mapply(MutualPartitionInfoSplits, rep(splits1, each=length(splits2)), splits2),
             length(splits2), length(splits1), dimnames = list(names(tree2), names(tree1)))
    }
  }
}

#' @describeIn MutualArborealInfo Takes splits instead of trees
#' @param splits1,splits2 Splits [#TODO document properly]
#' @export
MutualArborealInfoSplits <- function (splits1, splits2, reportMatching = FALSE) {
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  lnUnrootedN <- LnUnrooted.int(nTerminals)
  
  if (dimSplits1[2] < dimSplits2[2]) {
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
    
    tmp <- dimSplits1
    dimSplits1 <- dimSplits2
    dimSplits2 <- tmp
  }
  
  splits2 <- unname(splits2[rownames(splits1), , drop=FALSE])
  splits1 <- unname(splits1) # split1[split2] faster without names
  
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  inSplit1 <- colSums(splits1)
  inSplit2 <- colSums(splits2)
  notInSplit1 <- nTerminals - inSplit1 # TODO delete, and remove where used below?
  notInSplit2 <- nTerminals - inSplit2
  
  OneOverlap <- function(A1, A2) {
    if (A1 == A2) {
      # Return:
      LnRooted.int(A1) + LnRooted.int(nTerminals - A2)
    } else {
      if (A1 < A2) {
        tmp <- A2
        A2 <- A1
        A1 <- tmp
      }
      # Return:
      LnRooted.int(A1) + LnRooted.int(nTerminals - A2) - LnRooted.int(A1 - A2 + 1L) 
    }
  }
  
  pairScores <- matrix((mapply(function(i, j) {
    split1 <- splits1[, i]
    split2 <- splits2[, j]
    
    if (all(oneAndTwo <- split1[split2]) ||
        all(notOneNotTwo <- !split1[!split2])) {
      OneOverlap(inSplit1[i], inSplit2[j])
      
    } else if (all(notOneAndTwo <- !split1[split2]) ||
               all(oneNotTwo <- split1[!split2])) {
      OneOverlap(inSplit1[i], notInSplit2[j])
      
    } else {
      #in1 <- inSplit1[i]
      #out1 <- notInSplit1[i]
      #in2 <- inSplit2[j]
      #out2 <- notInSplit2[j]
      #splitTwoMoreEven <- min(in1, out1) < min(in2, out2)
      #contradictions <- c(min(sum(oneAndTwo), sum(notOneNotTwo)), 
      #                    min(sum(notOneAndTwo), sum(oneNotTwo)))
      #contradictionToFix <- which.min(contradictions)
      #
      #infoGainedBySwap <- if (contradictionToFix == 1L) {
      #  OneOverlap(in1, out2)
      #} else {
      #  OneOverlap(in1, in2)
      #}
      #
      #contradictionSize <- contradictions[contradictionToFix]
      #infoLostBySwap <- if (splitTwoMoreEven) {
      #  lchoose(in1, contradictionSize) + lchoose(out1, contradictionSize)
      #} else {
      #  lchoose(in2, contradictionSize) + lchoose(out2, contradictionSize)
      #}
      #c(SwapGain = infoGainedBySwap, SwapLoss = infoLostBySwap,
      #  SwapNet = infoGainedBySwap - infoLostBySwap)
      #
      #
      #contradictions <- c(sum(oneAndTwo), sum(notOneNotTwo), 
      #                    sum(notOneAndTwo), sum(oneNotTwo))
      lnUnrootedN
    }
  }, rep(seq_len(nSplits1), each=nSplits2), seq_len(nSplits2)
  ) - lnUnrootedN) / -log(2), nSplits2, nSplits1)
  
  if (nSplits2 == 1) {
    max(pairScores)
  } else {
    optimalMatching <- solve_LSAP(pairScores, TRUE)
    
    # Return:
    ret <- sum(pairScores[matrix(c(seq_along(optimalMatching), optimalMatching), ncol=2L)])
    if (reportMatching) {
      attr(ret, 'matching') <- optimalMatching
      attr(ret, 'pairScores') <- pairScores
      ret
    } else {
      ret
    }
  }
}

#' @describeIn VariationOfInfo Takes splits instead of trees
#' @param partitionQualityIndex Output of [SplitPairingInformationIndex] for
#' `n` taxa; calculated automatically if not specified, but passing a cached
#' value may improve performance.
#' @param bestMatchOnly Logical specifying whether to return how informative
#'  each split is about its best match only (`TRUE`) 
#'  or how informative each split is about each other split (`FALSE`).
#' @template splits12params
#' @export
VariationOfSplitPartitionInfo <- function (
  # TODO RENAME this function: it's a symmetric value saying how informative 
  # the partitions in one tree are about the best-matching partition in the 
  # other.
    splits1, splits2, reportMatching = FALSE,
    bestMatchOnly = TRUE,
    partitionQualityIndex = SplitPairingInformationIndex(dim(splits1)[1])
  ) {
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  lnUnrootedN <- LnUnrooted.int(nTerminals)
  
  if (dimSplits1[2] < dimSplits2[2]) {
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
    
    tmp <- dimSplits1
    dimSplits1 <- dimSplits2
    dimSplits2 <- tmp
  }
  
  splits2 <- unname(splits2[rownames(splits1), , drop=FALSE])
  splits1 <- unname(splits1) # split1[split2] faster without names
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  inSplit1 <- colSums(splits1)
  inSplit2 <- colSums(splits2)
  notInSplit1 <- nTerminals - inSplit1
  notInSplit2 <- nTerminals - inSplit2
  
  logTrees1 <- LnRooted(inSplit1) + LnRooted(notInSplit1)
  logTrees2 <- LnRooted(inSplit2) + LnRooted(notInSplit2)
  
  pairScores <- matrix((mapply(function(i, j) {
    split1 <- splits1[, i]
    split2 <- splits2[, j]
    
    variationOfInfo <- SplitEntropy(split1, split2)['vI']
    partitionQualityIndex[as.character(round(variationOfInfo, 6L))]
  }, rep(seq_len(nSplits1), each=nSplits2), seq_len(nSplits2)
  )), nSplits2, nSplits1)
  
  if (bestMatchOnly) {
    if (nSplits2 == 1) {
      min(pairScores)
    } else {
      optimalMatching <- solve_LSAP(pairScores, TRUE)
      
      # Return:
      ret <- sum(pairScores[matrix(c(seq_along(optimalMatching), optimalMatching), ncol=2L)])
      if (reportMatching) {
        attr(ret, 'matching') <- optimalMatching
        attr(ret, 'pairScores') <- pairScores
        ret
      } else {
        ret
      }
    }
  } else {
    # Return:
    sum(pairScores)
  }
}

#' @describeIn VariationOfInfo Takes splits instead of trees
#' @template splits12params
#' @export
VariationOfSplitArborealInfo <- function (splits1, splits2, reportMatching = FALSE) {
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  lnUnrootedN <- LnUnrooted.int(nTerminals)
  
  if (dimSplits1[2] < dimSplits2[2]) {
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
    
    tmp <- dimSplits1
    dimSplits1 <- dimSplits2
    dimSplits2 <- tmp
  }
  
  splits2 <- unname(splits2[rownames(splits1), , drop=FALSE])
  splits1 <- unname(splits1) # split1[split2] faster without names
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  inSplit1 <- colSums(splits1)
  inSplit2 <- colSums(splits2)
  notInSplit1 <- nTerminals - inSplit1
  notInSplit2 <- nTerminals - inSplit2
  
  logTrees1 <- LnRooted(inSplit1) + LnRooted(notInSplit1)
  logTrees2 <- LnRooted(inSplit2) + LnRooted(notInSplit2)
  
  OneOverlap <- function(A1, A2) {
    if (A1 == A2) {
      # Return:
      LnRooted.int(A1) + LnRooted.int(nTerminals - A2)
    } else {
      if (A1 < A2) {
        tmp <- A2
        A2 <- A1
        A1 <- tmp
      }
      # Return:
      LnRooted.int(A1) + LnRooted.int(nTerminals - A2) - LnRooted.int(A1 - A2 + 1L) 
    }
  }
  
  pairScores <- matrix((mapply(function(i, j) {
    split1 <- splits1[, i]
    split2 <- splits2[, j]
    
    logMutualTrees <- 
      if (all(oneAndTwo <- split1[split2]) ||
          all(notOneNotTwo <- !split1[!split2])) {
        OneOverlap(inSplit1[i], inSplit2[j])
        
      } else if (all(notOneAndTwo <- !split1[split2]) ||
                 all(oneNotTwo <- split1[!split2])) {
        OneOverlap(inSplit1[i], notInSplit2[j])
        
      } else {
        lnUnrootedN
      }
    logTrees1[i] + logTrees2[j] - logMutualTrees - logMutualTrees
  }, rep(seq_len(nSplits1), each=nSplits2), seq_len(nSplits2)
  )) / -log(2), nSplits2, nSplits1)
  
  
  
  if (nSplits2 == 1) {
    min(pairScores)
  } else {
    optimalMatching <- solve_LSAP(pairScores, FALSE)
    
    # Return:
    ret <- sum(pairScores[matrix(c(seq_along(optimalMatching), optimalMatching), ncol=2L)])
    if (reportMatching) {
      attr(ret, 'matching') <- optimalMatching
      attr(ret, 'pairScores') <- pairScores
      ret
    } else {
      ret
    }
  }
}

#' @describeIn MutualArborealInfo Takes splits instead of trees
#' @inheritParams MutualArborealInfoSplits
#' @export
MutualPartitionInfoSplits <- function (splits1, splits2, reportMatching = FALSE) {
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  
  if (dimSplits1[2] < dimSplits2[2]) {
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
    
    tmp <- dimSplits1
    dimSplits1 <- dimSplits2
    dimSplits2 <- tmp
  }
  
  splits2 <- unname(splits2[rownames(splits1), , drop=FALSE])
  splits1 <- unname(splits1) # split1[split2] faster without names
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  
  pairScores <- -log2(matrix(mapply(function(i, j) {
    split1 <- splits1[, i]
    split2 <- splits2[, j]
    
    SplitMatchProbability(split1, split2)
  }, rep(seq_len(nSplits1), each=nSplits2), seq_len(nSplits2))
  , nSplits2, nSplits1))
  
  if (nSplits2 == 1) {
    max(pairScores)
  } else {
    optimalMatching <- solve_LSAP(pairScores, TRUE)
    # Previously:
    # sum(pairScores[cbind(seq_along(optimalMatching), optimalMatching)])
    # Now (30% faster):
    ret <- sum(pairScores[matrix(c(seq_along(optimalMatching), optimalMatching), ncol=2L)])
    
    if (reportMatching) {
      attr(ret, 'matching') <- optimalMatching
      attr(ret, 'pairScores') <- pairScores
      ret
    } else {
      ret
    }
  }
}

#' Are splits compatible?
#' 
#'
#' Splits are compatible if they are concave; i.e. they can both be true
#' simultaneously.
#' @template split12params
#' 
#' @author Martin R. Smith
#' @export
SplitsCompatible <- function (split1, split2) {
  all (split1[split2]) ||
    all(split1[!split2]) ||
    all(!split1[split2]) ||
    all(!split1[!split2])
}
