#' Number of trees matching a bipartition split
#' 
#' Calculates the number of unrooted bifurcated trees that are consistent with
#' a bipartition split that divides taxa into groups of size `A` and `B`.
#' 
#' @param A,B Number of taxa in each partition.
#' 
#' @author Martin R. Smith
#' 
#' @family split information functions
#' @export
TreesMatchingSplit <- function (A, B) {
  if (A == 0) NUnrooted(B) else
  if (B == 0) NUnrooted(A) else
  NRooted(A) * NRooted(B)
}

#' @describeIn TreesMatchingSplit Logarithm of the number of trees matching a split.
#' @export
LogTreesMatchingSplit <- function (A, B) {
  if (A == 0) LnUnrooted.int(B) else
  if (B == 0) LnUnrooted.int(A) else
  LnRooted.int(A) + LnRooted.int(B)
}

#' Character information content
#' 
#' Calculates the phylogenetic information content of a given character.
#' 
#' @param tokens Character vector specifying the tokens assigned to each taxon for 
#' a character.  Example: `c(0, 0, 0, 1, 1, 1, '?', '-')`.
#' 
#' Note that ambiguous tokens such as `(01)` are not supported, and should be
#' replaced with `?`.`
#' 
#' @return Phylogenetic information content of the character, in bits.
#' 
#' @author Martin R. Smith
#' @export
CharacterInformation <- function (tokens) {
  tokenCounts <- table(tokens)
  # Our character splits our taxa into groups with the same token
  # ?s and -s are best ignored
  splits <- tokenCounts[!(names(tokenCounts) %in% c('?', '-'))]
  
  # Information content = -log2(probability)
  # Probability of a tree being consistent with our character is
  # n trees consistent with character / n trees with that many tips
  # NUnrootedMult(splits) / NUnrooted(splits)
  # As we are working with large numbers we can use logarithms
  lnP <- LnUnrootedMult(splits) - LnUnrooted(sum(splits))
  log2P <- lnP / log(2)
  information <- -log2P
  
  # Return: 
  information
}

#' Information content of a split
#' 
#' `SplitInformation` calculates the information content of a split, based on 
#' the entropy of the subset of trees consistent with the split; a split that
#' is consistent with a smaller number of trees will have a higher information
#' content.
#'
#' @inheritParams TreesMatchingSplit
#'
#' @return Information content of the split, in bits.
#' 
#' @examples 
#'   # Eight tips can be split evenly:
#'   SplitInformation (4, 4)
#'   
#'   # or unevenly, which is less informative:
#'   SplitInformation (2, 6)
#' 
#' @author Martin R. Smith
#' @family split information functions
#' @export
SplitInformation <- function (A, B) {
  -(LogTreesMatchingSplit(A, B) - LnUnrooted.int(A + B)) / log(2)
}

#' @describeIn SplitInformation Information content of a multi-partition split.
#' @param partitionSizes Integer vector specifying the number of taxa in each 
#' partition of a multi-partition split.
#' @export
MultiSplitInformation <- function (partitionSizes) {
  -(LnUnrootedMult(partitionSizes) - LnUnrooted.int(sum(partitionSizes))) / log(2)
}

#' Distributions of taxa consistent with a partition pair.
#' 
#' Number of terminal arrangements matching a specified configuration of 
#' two partitions.
#' 
#' Consider partitions that divide eight terminals, labelled A to H.
#' 
#' \tabular{rcll}{
#'   Bipartition 1:\tab ABCD:EFGH\tab A1 = ABCD\tab B1 = EFGH \cr
#'   Bipartition 2:\tab ABE:CDFGH\tab A2 = ABE\tab B2 = CDFGH
#' }
#' 
#' This can be represented by an association matrix:
#' 
#' \tabular{rll}{
#'      \tab *A2* \tab *B2* \cr
#' *A1* \tab AB  \tab C   \cr
#' *B1* \tab E   \tab FGH
#' }
#' 
#' The cells in this matrix contain 2, 1, 1 and 3 terminals respectively; this
#' four-element vector (`c(2, 1, 1, 3)`) is the `configuration` implied by 
#' this pair of bipartition splits.
#' 
#' @param configuration Integer vector of length four specifying the number of 
#' terminals that occur in both 
#'   (1) splits A1 and A2; 
#'   (2) splits A1 and B2;
#'   (3) splits B1 and A2;
#'   (4) splits B1 and B2.
#'   
#' @return The number of ways to distribute `sum(configuration)` taxa according
#'  to the specified pattern.
#'
#' @examples 
#' NPartitionPairs(c(2, 1, 1, 3))
#'  
#' @author Martin R. Smith
#' @export
NPartitionPairs <- function (configuration) {
  choose(sum(configuration[c(1, 3)]), configuration[1]) *
  choose(sum(configuration[c(2, 4)]), configuration[2])
}

#' @describeIn SplitMutualInformation Number of trees consistent with two splits.
#' @family split information functions
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

#' @describeIn SplitMutualInformation Natural logarithm of `TreesConsistentWithTwoSplits`.
#' @export
LogTreesConsistentWithTwoSplits <- function (n, A1, A2=A1) {
  smallSplit <- min(A1, A2)
  bigSplit <- max(A1, A2)
  
  if (smallSplit == 0) return (LogTreesMatchingSplit(bigSplit, n - bigSplit))
  if (bigSplit == n) return (LogTreesMatchingSplit(smallSplit, n - smallSplit))
  
  # Return:
  LnRooted(bigSplit - smallSplit + 1L) + LnRooted(smallSplit) + LnRooted(n - bigSplit)
}

#' Number of trees consistent with split
#' 
#' Calculates the number of unrooted bifurcating trees consistent with the 
#' specified multi-partition split, using the formula of Carter _et al_. (1990).
#' 
#' @template splitsParam
#' 
#' @return `UnrootedTreesMatchingSplit` returns an integer specifying the
#'  number of unrooted bifurcating trees consistent with the specified split.
#' 
#'
#' @examples 
#'  UnrootedTreesMatchingSplit(c(3, 5))
#'  UnrootedTreesMatchingSplit(c(3, 2, 1, 2))
#' 
#' @references 
#' \insertRef{Carter1990}{TreeSearch}, Theorem 2.
#'
#' @author Martin R. Smith
#' @family split information functions
#' @export
UnrootedTreesMatchingSplit <- function (splits) {
  splits <- splits[splits > 0L]
  totalTips <- sum(splits)
  tipsMinusLengthSplits <- totalTips - length(splits)
  # use exp and log as it's just as fast, but less likely to overflow to Inf
  exp(sum(LogDoubleFactorial(totalTips + totalTips - 5L),
          LogDoubleFactorial(splits + splits - 3L)) -
        LogDoubleFactorial(tipsMinusLengthSplits + tipsMinusLengthSplits - 1L))
}
