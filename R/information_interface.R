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
  if (A == 0) LnUnrooted(B) else
  if (B == 0) LnUnrooted(A) else
  LnRooted(A) + LnRooted(B)
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
  (LogTreesMatchingSplit(A1, n - A1) + LogTreesMatchingSplit(A2, n - A2) -
    LnUnrooted(n) - LogTreesConsistentWithTwoSplits(n, A1, A2)) / -log(2)
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
