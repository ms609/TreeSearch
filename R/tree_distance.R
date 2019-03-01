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
#' @references 
#'  \insertRef{SmithDist}{TreeSearch}
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

#' Variation of information between two trees
#' 
#' See Melia 2007... #TODO document
#'
#' @inheritParams MutualArborealInfo
#' 
#' 
#' 
#' @references {
#'   \insertRef{Meila2007}{TreeSearch}
#' }
#'
#' @author Martin R. Smith
#' @export
#' Tree distance based on joint information content of splits
#' 
#' #TODO Needs documenting and describing fully.
#' 
#' @inheritParams MutualArborealInfo
#' @inheritSection references MutualArborealInfo 
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
  
  if (dim(splits1)[2] < dim(splits2)[2]) {
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
  }
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  lnUnrootedN <- LnUnrooted.int(nTerminals)
  
  splits2 <- unname(splits2[rownames(splits1), , drop=FALSE])
  splits1 <- unname(splits1) # split1[split2] faster without names
  
  
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  
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
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  inSplit1 <- colSums(splits1)
  inSplit2 <- colSums(splits2)
  notInSplit2 <- nTerminals - inSplit2
  
  pairScores <- matrix((mapply(function(i, j) {
    split1 <- splits1[, i]
    split2 <- splits2[, j]
    
    if (all(split1[split2]) || all(!split1[!split2])) {
      OneOverlap(inSplit1[i], inSplit2[j])
      
    } else if (all(!split1[split2]) || all(split1[!split2])) {
      OneOverlap(inSplit1[i], notInSplit2[j])
      
    } else {
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
#' @template splits12params
#' @export
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
