#' Number of trees matching a bipartition split
#' 
#' @param A,B Number of taxa on each side of split
#' 
#' @author Martin R. Smith
#' 
#' @concept Split information
#' @export
TreesMatchingSplit <- function (A, B) {
  if (A == 0) NUnrooted(B) else
  if (B == 0) NUnrooted(A) else
  NRooted(A) * NRooted(B)
}

#' @describeIn TreesMatchingSplit Logarithm of Trees Matching Split
#' @export
LogTreesMatchingSplit <- function (A, B) {
  if (A == 0) LnUnrooted.int(B) else
  if (B == 0) LnUnrooted.int(A) else
  LnRooted.int(A) + LnRooted.int(B)
}

#' Mutual information of two splits
#' 
#' @param n Integer specifying the number of terminals.
#' @param A1,A2 Integers specifying the number of taxa in A1 and A2, 
#' once the splits have been arranged such that all taxa in A1 are also
#' in A2, or all taxa in A2 are also in A1.
#' 
#' @return The information that two splits have in common, in bits.
#' 
#' @author Martin R. Smith
#' @concept Split information
#' @export
MutualInformation <- function(n, A1, A2=A1) {
  (LogTreesMatchingSplit(A1, n - A1) 
   + LogTreesMatchingSplit(A2, n - A2)
   - LogTreesConsistentWithTwoSplits(n, A1, A2) 
   - LnUnrooted(n)) / -log(2)
}

#' Information content of a split
#'
#' @inheritParams TreesMatchingSplit
#'
#' @return Information content of the split in bits.
#' @author Martin R. Smith
#' @concept Split information
#' @export
SplitInformation <- function (A, B) {
  -log2(TreesMatchingSplit(A, B) / NUnrooted(A+B))
}

#' @describeIn SplitInformation Information content of a multi-partition split.
#' @param partitionSizes Integer vector specifying the number of taxa in each 
#' partition of a multi-partition split.
#' @export
MultiSplitInformation <- function (partitionSizes) {
  -log2(NUnrootedMult(partitionSizes) / NUnrooted(sum(partitionSizes)))
}


#' Probability that two random splits will be at least as similar as these two
#' @template split1Param
#' @template split2Param
#' 
#' @return The natural logarithm of the probability of observing two splits,
#' spliting the terminals into bipartitions of the sizes given,
#'  that match as well as `split1` and `split2` do.
LnSplitMatchProbability <- function (split1, split2) {
  
  if (length(split1) != length(split2)) stop("Split lengths differ")
  
  A1A2 <- sum(split1 & split2)
  A1B2 <- sum(split1 & !split2)
  B1A2 <- sum(!split1 & split2)
  B1B2 <- sum(!split1 & !split2)
  c(A1A2 = A1A2, A1B2 = A1B2, B1A2 = B1A2, B1B2 = B1B2)
  
  A1 <- A1A2 + A1B2
  A2 <- A1A2 + B1A2
  B1 <- B1A2 + B1B2
  B2 <- A1B2 + B1B2
  n <- A1 + B1
  paste0(A1, ':', B1, "; ", A2, ':', B2)
  
  # Return:
  if (A1A2 == 0 || B1B2 == 0) {
    # Overlapping groups = A1-B2, B1-A2
    if (A1 > B2) {
      lchoose(A1, B2) - lchoose(n, B2)
    } else {
      lchoose(B2, A1) - lchoose(n, A1)
    }
  } else if (A1B2 == 0 || B1A2 == 0) {
    # Overlapping groups = A1-A2, B1-B2
    if (A1 > A2) {
      lchoose(A1, A2) - lchoose(n, A2)
    } else {
      lchoose(A2, A1) - lchoose(n, A1)
    }
  } else {
    Ways <- function (bigMatch, smallMatch, matchOverlap) {
      chosenFromBigMatch <- matchOverlap:smallMatch
      log(sum(choose(bigMatch, chosenFromBigMatch) * 
                choose(n - bigMatch, smallMatch - chosenFromBigMatch))) -
            lchoose(n, smallMatch)
    }
    
    # Clades discordant
    if (A1A2 > A1B2 || B1B2 > B1A2) {
      Ways(max(A1, A2), min(A1, A2), A1A2)
    # == Ways(max(B1, B2), min(B1, B2), B1B2)
    } else {
      Ways(max(A1, B2), min(A1, B2), A1B2)
    # == Ways(max(B1, A2), min(B1, A2), B1A2)
    }

  }
}

SplitEntropy <- function (split1, split2=split1) {
  A1A2 <- sum(split1 & split2)
  A1B2 <- sum(split1 & !split2)
  B1A2 <- sum(!split1 & split2)
  B1B2 <- sum(!split1 & !split2)
  overlaps <- c(A1A2, A1B2, B1A2, B1B2)
  
  
  A1 <- A1A2 + A1B2
  A2 <- A1A2 + B1A2
  B1 <- B1A2 + B1B2
  B2 <- A1B2 + B1B2
  n <- A1 + B1
  
  H <- function(p) -sum(p*log(p))
  h1 <- H(c(A1, B1) / n)
  h2 <- H(c(A2, B2) / n)
  jointH <- H(overlaps[overlaps > 0L] / n)
  mutualInformation <- h1 + h2 - jointH
  
  # Return:
  c(h1 = h1, h2 = h2, jointH = jointH, i = mutualInformation)  
}

#' Joint Information of two splits
#'
#' Because some information is common to both splits (`MutualInformation`),
#' the joint information of two splits will be less than the sum of the
#' information of the splits taken separately -- unless the splits are
#' contradictory.
#' 
#' Split Y1 is defined as dividing taxa into the two sets A1 and B1,
#' and Y2=A2:B2.
#'
#' @param A1A2,A1B2,B1A2,B1B2 Number of taxa in splits A1 and A2 (etc.)
#' @author Martin R. Smith
#' @concept Split information
#' @export
JointInformation <- function(A1A2, A1B2, B1A2, B1B2) {
  A1 <- A1A2 + A1B2
  B1 <- B1A2 + B1B2
  A2 <- A1A2 + B1A2
  B2 <- A1B2 + B1B2
  n  <- A1 + B1
  
  mutualInformation <- if(any(c(A1A2, A1B2, B1A2, B1B2) == 0)) {
    -log2(TreesConsistentWithTwoSplits(n, A1, A2) / NUnrooted(n))
  } else {
    0
  }
  
  # Return:
  SplitInformation(A1, B1) + SplitInformation(A2, B2) - mutualInformation
}

#' Restricted joint information
#' 
#' #TODO keep and describe, or delete.
#' 
#' Information based on the proportion of trees that are consistent with
#' Y1 OR Y2, AND consistent with the information that Y1 and Y2 have in common.
#'  
#' @inheritParams JointInformation
#' @author Martin R. Smith
#' @concept Split information
#' @export
FullMutualInformation <- function(A1A2, A1B2, B1A2, B1B2) {
  A1 <- A1A2 + A1B2
  B1 <- B1A2 + B1B2
  A2 <- A1A2 + B1A2
  B2 <- A1B2 + B1B2
  n  <- A1 + B1
  
  
  BuildSplitFromCommonInformation <- function (A1, A1a, B1, B1a) {
    A1b <- A1 - A1a
    B1b <- B1 - B1a
    
    lnSubtreesA <- if (A1a == 0 || A1b == 0) {
      LnUnrooted(A1) 
    } else {
      LnRooted(A1a) + LnRooted(A1b)
    }
    
    lnSubtreesB <- if (B1a == 0 || B1b == 0) {
      LnUnrooted(B1) 
    } else {
      LnRooted(B1a) + LnRooted(B1b)
    }
    
    lnTrees <- lnSubtreesA + log(A1 + A1 - 3L) + 
      lnSubtreesB + log(B1 + B1 - 3L) - 
      LnUnrooted(A1 + B1)
    
    # Return: 
    lnTrees / -log(2)
  }
  
  
  # Return:
  if(any(c(A1A2, A1B2, B1A2, B1B2) == 0)) {
    a1Pair <- if (A1A2 == 0 || B1B2 == 0) B2 else A2
    # Splits are consistent
    MutualInformation(n, A1, a1Pair)
  } else {
    # Splits are inconsistent
    SplitInformation(A1A2, B1B2) +
      SplitInformation(B1A2, A1B2)
  }
  
}

#' Number of trees consistent with two splits
#'
#' @inheritParams MutualInformation
#' @author Martin R. Smith
#' @concept Split information
#' @export
TreesConsistentWithTwoSplits <- function (n, A1, A2=A1) {
  
  smallSplit <- min(A1, A2)
  bigSplit <- max(A1, A2)
  
  if (smallSplit == 0) return (TreesMatchingSplit(bigSplit, n - bigSplit))
  if (bigSplit == n) return (TreesMatchingSplit(smallSplit, n - smallSplit))
  
  overlap <- bigSplit - smallSplit
  
  #  Here are two spits:
  #  AA OOO BBBBBB
  #  11 111 000000
  #  11 000 000000
  #  
  #  There are (2O - 5)!! unrooted trees of the overlapping taxa (O)
  #  There are (2A - 5)!! unrooted trees of A
  #  There are (2B - 5)!! unrooted trees of B
  #  
  #  There are 2A - 3 places on A that A can be attached to O.
  #  There are 2O - 3 places on O to which A can be attached.
  #  
  #  There are 2B - 3 places in B that B can be attached to (O+A).
  #  There are 2O - 3 + 2 places on O + A to which B can be attached:
  #   2O - 3 places on O, plus the two new edges created when A was joined to O.
  #  
  #  (2A - 3)(2A - 5)!! == (2A - 3)!!
  #  (2B - 3)(2B - 5)!! == (2B - 3)!!
  #  (2O - 1)(2O - 3)(2O - 5)!! == (2O - 1)!!
  #  
  #  We therefore want NRooted(A) * NRooted(O + 1) * NRooted(B)
  #  
  #  O = overlap = bigSplit - smallSplit
  #  bigSplit - overlap = smallSplit = either A or B
  #  n - bigSplit = either B or A
  
  
  # Return:
  NRooted(overlap + 1L) * 
    NRooted(smallSplit) *
    NRooted(n - bigSplit)
}

#' @describeIn TreesConsistentWithTwoSplits Logarithm thereof
#' @export
LogTreesConsistentWithTwoSplits <- function (n, A1, A2=A1) {
  smallSplit <- min(A1, A2)
  bigSplit <- max(A1, A2)
  
  if (smallSplit == 0) return (LogTreesMatchingSplit(bigSplit, n - bigSplit))
  if (bigSplit == n) return (LogTreesMatchingSplit(smallSplit, n - smallSplit))
  
  # Return:
  LnRooted(bigSplit - smallSplit + 1L) + LnRooted(smallSplit) + LnRooted(n - bigSplit)
}

#' Number of unrooted trees consistent with splits
#' @template splitsParam
#' 
#' @references 
#' \insertRef{Carter1990}{TreeSearch}, Theorem 2.
#'
UnrootedTreesMatchingSplit <- function (splits) {
  splits <- splits[splits > 0L]
  totalTips <- sum(splits)
  tipsMinusLengthSplits <- totalTips - length(splits)
  # use exp and log as it's just as fast, but less likely to overflow to Inf
  exp(sum(LogDoubleFactorial(totalTips + totalTips - 5L),
          LogDoubleFactorial(splits + splits - 3L)) -
        LogDoubleFactorial(tipsMinusLengthSplits + tipsMinusLengthSplits - 1L))
}
