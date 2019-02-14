#' Information-based generalized Robinson-Fould distance between two trees
#'
#' @param tree1,tree2 Trees of class `phylo`, with tips labelled identically.
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
  if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
    stop("Tree tips must bear identical labels")
  }

  if (tree1$Nnode < tree2$Nnode) {
    splits1 <- Tree2Splits(tree2)
    splits2 <- Tree2Splits(tree1)
  } else {
    splits1 <- Tree2Splits(tree1)
    splits2 <- Tree2Splits(tree2)
  }

  nTerminals <- nrow(splits1)
  lnUnrootedN <- LnUnrooted(nTerminals)

  splits2 <- splits2[rownames(splits1), , drop=FALSE]

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
  
  
  pairScores <- (apply(splits1, 2, function (split1) {
    apply(splits2, 2, function (split2) {
      if (all(split1[split2]) || all(!split1[!split2])) {
        OneOverlap(sum(split1), sum(split2))
        
      } else if (all(!split1[split2]) || all(split1[!split2])) {
        OneOverlap(sum(split1), sum(!split2))
        
      } else {
        lnUnrootedN
      }
    })
  }) - lnUnrootedN) / -log(2)

  if (is.null(dim(pairScores))) {
    # Only one split in splits2, so apply returns a vector instead of an array
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
