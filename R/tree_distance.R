#' Information-based generalized Robinson-Fould distance between two trees
#'
#' @param tree1,tree2 Trees of class `phylo`, with tips labelled identically,
#' or lists of such trees to undergo pairwise comparison.
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
InfoTreeDist <- function (tree1, tree2) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      MutualSplitInformation(Tree2Splits(tree1), Tree2Splits(tree2))
    } else {
      splits1 <- Tree2Splits(tree1)
      vapply(tree2, 
             function (tr2) MutualSplitInformation(splits1, Tree2Splits(tr2)),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- Tree2Splits(tree2)
      vapply(tree1, 
             function (tr2) MutualSplitInformation(splits1, Tree2Splits(tr2)),
             double(1))
    } else {
      splits1 <- lapply(tree1, Tree2Splits)
      splits2 <- lapply(tree2, Tree2Splits)
      matrix(mapply(MutualSplitInformation, rep(splits1, each=length(splits2)), splits2),
             length(splits2), length(splits1), dimnames = list(names(tree2), names(tree1)))
    }
  }
}

#' @describeIn InfoTreeDist Takes splits instead of trees
#' @param splits1,splits2 Splits [#TODO document properly]
#' @export
MutualSplitInformation <- function (splits1, splits2) {
  
  if (dim(splits1)[2] < dim(splits2)[2]) {
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
  }
  
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  lnUnrootedN <- LnUnrooted(nTerminals)
  
  splits2 <- splits2[rownames(splits1), , drop=FALSE]
  
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  
  OneOverlap <- function(A1, A2) {
    if (A1 == A2) {
      # Return:
      LnRooted(A1) + LnRooted(nTerminals - A2)
    } else {
      if (A1 < A2) {
        tmp <- A2
        A2 <- A1
        A1 <- tmp
      }
      # Return:
      LnRooted(A1) + LnRooted(nTerminals - A2) - LnRooted(A1 - A2 + 1L) 
    }
  }
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  pairScores <- matrix((mapply(function(i, j) {
    split1 <- splits1[, i]
    split2 <- splits2[, j]
    if (all(split1[split2]) || all(!split1[!split2])) {
      OneOverlap(sum(split1), sum(split2))
      
    } else if (all(!split1[split2]) || all(split1[!split2])) {
      OneOverlap(sum(split1), sum(!split2))
      
    } else {
      lnUnrootedN
    }
  }, rep(seq_len(nSplits1), each=nSplits2), seq_len(nSplits2)
  ) - lnUnrootedN) / -log(2), nSplits2, nSplits1)
  
  if (nSplits2 == 1) {
    max(pairScores)
  } else {
    optimalMatching <- solve_LSAP(pairScores, TRUE)
    sum(pairScores[cbind(seq_along(optimalMatching), optimalMatching)])
  }
  
}

#' Are splits compatible?
#' 
#' @param split1,split2 Logical vectors listing terminals in same order, such that
#' each terminal is identified as a member of the ingroup (`TRUE`) or outgroup 
#' (`FALSE`) of a bipartition split.
#'
#' Splits are compatible if they are concave; i.e. they can both be true
#' simultaneously.
SplitsCompatible <- function (split1, split2) {
  all (split1[split2]) ||
    all(split1[!split2]) ||
    all(!split1[split2]) ||
    all(!split1[!split2])
}
