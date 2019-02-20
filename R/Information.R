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
#' @param n Number of terminals
#' @param A1,A2 Number of terminals on overlapping side of each split.
#' 
#' @return The information that two splits have in common
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
#' @return Information content in bits.
#' @author Martin R. Smith
#' @concept Split information
#' @export
SplitInformation <- function (A, B) {
  -log2(TreesMatchingSplit(A, B) / NUnrooted(A+B))
}

#' Probability that two random splits will be at least as similar as these two
#' @template split1Param
#' @template split2Param
#' 
#' @return The natural logarithm of the probability of observing two splits,
#' of the sizes input, that match as well as `split1` and `split2` do.
LnSplitMatchProbability <- function (split1, split2) {
  
  # Define side A as the side that contains the first taxon
  if (!split1[1]) split1 <- !split1
  if (!split2[1]) split2 <- !split2
  
  A1A2 <- sum(split1 & split2) - 1L # -1 degree of freedom because A is in both
  A1B2 <- sum(split1 & !split2)
  B1A2 <- sum(!split1 & split2)
  B1B2 <- sum(!split1 & !split2)
  
  A1 <- A1A2 + A1B2
  A2 <- A1A2 + B1A2
  
  # Return:
  if (A2 > A1) {
    B2 <- A1B2 + B1B2
    log(sum(choose(A2, (A1A2:A1)) * 
              choose(B2, A1 - (A1A2:A1)))) - lchoose(n - 1L, A1)
  } else {
    B1 <- B1A2 + B1B2
    log(sum(choose(A1, (A1A2:A2)) * 
              choose(B1, A2 - (A1A2:A2)))) - lchoose(n - 1L, A2)
  }
  ## Return:
  #if (A2 > A1) {
  #  lchoose(A2, A1) - lchoose(n, A1)
  #} else {
  #  lchoose(A1, A2) - lchoose(n, A2)
  #}
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
