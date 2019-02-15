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

#' Joint Information of two splits
#'
#' Because some information is common to both splits (`MutualInformation`),
#' the joint information of two splits will be less than the sum of the
#' information of the splits taken separately -- unless the splits are
#' contradictory.
#'
#' @inheritParams MutualInformation
#' @author Martin R. Smith
#' @concept Split information
#' @export
JointInformation <- function(n, A1, A2=A1) {
  -log2(TreesConsistentWithTwoSplits(n, A1, A2) / NUnrooted(n))
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
