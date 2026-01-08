#' Maximum likelihood estimate of consistency
#'
#' Compute MLCI for all columns of a character matrix.
#'
#' @param tree an ape::phylo tree
#' @param charmat matrix with rows = tips, columns = characters
#' @param precision stopping criterion on SE(MLCI)
#' @param maxResample maximum number of resamples
#' @param lower,upper bounds for t_ML
#' @param tol tolerance for t_ML optimizer
#'
#' @return data.frame with columns:
#'   t_hat, score, rmMean, rmSE, bestScore, mlci, nResample
#' @export
MLCI <- function(tree, charmat, precision = 1e-2, maxResample = 10000) {
  requireNamespace("TreeTools")
  
  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
  edge <- tree$edge
  Ntip <- length(tree$tip.label)
  
  # Compute best split once
  splits <- TreeTools::as.Splits(tree)
  splitImb <- TreeTools::SplitImbalance(splits)
  bestSplitIdx <- which.min(splitImb)
  bestSplit <- as.integer(as.logical(splits[[bestSplitIdx]]))
  
  # Helper: map tip states of a character column to 0..k-1 integer
  prep_tip_states_single <- function(col) {
    col <- as.character(col)
    lvl <- sort(unique(col))
    k <- length(lvl)
    map <- setNames(seq_len(k)-1L, lvl)
    tip_states <- matrix(as.integer(map[col]), nrow = Ntip, ncol = 1)
    list(tip_states = tip_states, k = k)
  }
  
  # Call C++ engine
  res <- MLCI_rcpp(
    edge = edge,
    tipStates = charmat,
    bestSplit = bestSplit,
    precision = precision,
    maxResample = maxResample
  )
  
  res
}
