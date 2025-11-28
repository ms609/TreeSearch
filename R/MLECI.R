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
MLCI <- function(tree, charmat, precision = 1e-2, maxResample = 10000, lower = 1e-8, upper = 10, tol = 1e-6, maxK = 10) {
  requireNamespace("TreeTools")
  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
  edge <- tree$edge
  Ntip <- length(tree$tip.label)
  nChars <- ncol(charmat)
  
  # prepare engine once
  eng <- mlci_make_engine(edge, Ntip, maxK)
  
  # compute best split once
  splits <- TreeTools::as.Splits(tree)
  im <- TreeTools::SplitImbalance(splits)
  best_split <- splits[[which.min(im)]]
  best_states <- as.integer(as.logical(best_split)) # 0/1
  
  # result container
  res <- data.frame(t_hat = numeric(nChars), score = numeric(nChars), rmMean = numeric(nChars),
                    rmSE = numeric(nChars), bestScore = numeric(nChars), mlci = numeric(nChars),
                    mlciSE = numeric(nChars), nResample = integer(nChars))
  
  for (j in seq_len(nChars)) {
    col <- as.character(charmat[, j])
    lvl <- sort(unique(col))
    k <- length(lvl)
    map <- setNames(seq_len(k) - 1L, lvl)
    obs_states <- as.integer(map[col])
    
    out <- mlci_resample_engine(eng, obs_states, best_states, precision, maxResample, lower, upper, tol)
    
    res$t_hat[j] <- out$t_hat_obs
    res$score[j] <- out$score
    res$rmMean[j] <- out$rmMean
    res$rmSE[j] <- out$rmSE
    res$bestScore[j] <- out$bestScore
    res$mlci[j] <- out$mlci
    res$mlciSE[j] <- out$mlciSE
    res$nResample[j] <- out$nResample
  }
  
  res
}
