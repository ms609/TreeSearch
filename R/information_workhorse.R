# Interfaces to this function are in information_interface.R

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


#' Double Factorial
#' 
#' @param n Vector of integers.
#' 
#' @return Returns the double factorial, n x (n - 2) x (n - 4) x (n - 6) x ...
#' 
#' @examples {
#' DoubleFactorial (0) # Return 1 if n < 2
#' DoubleFactorial (2) # 2
#' DoubleFactorial (5) # 1 x 3 x 5
#' DoubleFactorial (8) # 2 x 4 x 6 x 8
#' }
#' 
#' @author Martin R. Smith
#' @export
DoubleFactorial <- function (x) {
  if (any(x > 300)) stop("301!! is too large to store as an integer. Use LogDoubleFactorial instead.")
  
  x[x < 2] <- 1
  doubleFactorial[x]
  
  #
  #odds <- as.logical(x %% 2)
  #
  #oddX <- x[odds]
  #xPlusOneOverTwo <- (oddX + 1) / 2
  #evenX <- x[!odds]
  #xOverTwo <- evenX / 2
  #
  #ret <- integer(length(x))
  #ret[odds] <- gamma(oddX + 1L) / (gamma(xPlusOneOverTwo) * 2^(xPlusOneOverTwo - 1L))
  #ret[!odds] <- evenX * gamma(xOverTwo) * 2^(xOverTwo - 1L)
  #
  ## Return:
  #ret
}

# Memoizing this function makes it MUCH slower...
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
LogDoubleFactorial <- (function (x) {
  x[x < 2] <- 1
  if (all(x < 50000L)) {
    # Return from cache
    logDoubleFactorial[x]
  } else {
    
    odds <- as.logical(x %% 2)
    
    oddX <- x[odds]
    xPlusOneOverTwo <- (oddX + 1) / 2
    evenX <- x[!odds]
    xOverTwo <- evenX / 2
    
    ret <- integer(length(x))
    ret[odds] <- lgamma(oddX + 1L) - (lgamma(xPlusOneOverTwo) + (xPlusOneOverTwo - 1) * log(2))
    ret[!odds] <- log(evenX) + lgamma(xOverTwo) + (xOverTwo - 1) * log(2)
    
    # Return:
    ret
  }
})
