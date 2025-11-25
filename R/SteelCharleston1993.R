#' Number of tree labellings with _k_ steps
#' 
#' Number of ways to assign labels `0` and `1` to leaves of a rooted binary tree
#' such that the tree has Fitch parsimony length `k`.
#' As only two labels are available, this is independent of tree shape.
#' Theorem 1 in \insertCite{Steel1993}{TreeSearch}
#' 
#' @param tips Integer giving number of leaves; or tree whose leaves will be
#' counted with [NTip()].
#' @param steps Integer giving number of parsimony steps.
#' @param log Logical or numeric; if `TRUE`, results will be log transformed
#' in base $e$; if numeric, in base `log`.
#' @returns The number of possible binary leaf labellings.
#' @references \insertAllCited{}
#' @importFrom TreeTools NTip
#' @export
BinaryLabellings <- function(tips, steps, log = FALSE) {
  
  if (!is.numeric(tips)) {
    tips <- NTip(tips)
  }
  if (!is.numeric(tips)) {
    warning("`tips` must be numeric")
    return(0)
  }
  
  numerator <- tips - steps - 1
  if (isFALSE(log)) {
    twoK <- 2 ^ steps
    n0 <- choose(numerator, steps)
    n01 <- choose(numerator, steps - 1)
    # Return:
    twoK * sum(n0, n0, n01)
  } else {
    ln0 <- lchoose(numerator, steps)
    ln01 <- lchoose(numerator, steps - 1)
    lnX <- log(2) * steps + LogSumExp(ln0 + log(2), ln01)
    # Return:
    if (is.numeric(log)) {
      lnX / log(log)
    } else {
      lnX
    }
  }
}
