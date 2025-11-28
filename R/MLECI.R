#' Maximum likelihood estimate of consistency
#' 
#' Compute MLCI for a matrix of characters
#' @param tree an ape::phylo tree
#' @param charmat matrix-like with columns = characters, rows = tips (tokens)
#' @param precision desired standard error for the mean of rmScore (stop criterion)
#' @param maxResample maximum number of resamples to attempt
#' @param lower,upper search bounds for t in MLE
#' @param tol optimizer tolerance
#' @return data.frame with columns: t_hat, score, rmMean, rmSE, bestScore, mlci
#' @examples
#' tree <- TreeTools::BalancedTree(8)
#' chars <- cbind(c(0, 0, 0, 0, 1, 1, 1, 1),
#'                c(0, 0, 0, 0, 1, 0, 0, 0),
#'                c(0, 0, 0, 0, 0, 1, 1, 1),
#'                c(0, 1, 0, 1, 0, 1, 0, 1))
#'  MLCI(tree, chars)
#' 
#' @export
MLCI <- function(tree, charmat, precision = 1e-3, maxResample = 10000, lower = 1e-8, upper = 10, tol = 1e-6) {
  requireNamespace("TreeTools")
  if (!inherits(tree, "phylo")) stop("tree must be of class 'phylo'")
  edge <- tree$edge
  Ntip <- length(tree$tip.label)
  nEdge <- nrow(edge)
  nChars <- ncol(charmat)
  
  # Precompute splits (for bestScore) once
  splits <- TreeTools::as.Splits(tree)
  splitImb <- TreeTools::SplitImbalance(splits)
  bestSplitIdx <- which.min(splitImb)
  bestSplit <- splits[[bestSplitIdx]] # logical vector length Ntip
  
  # result containers
  res <- data.frame(t_hat = numeric(nChars), score = numeric(nChars), rmMean = numeric(nChars), rmSE = numeric(nChars), bestScore = numeric(nChars), mlci = numeric(nChars))
  
  for (j in seq_len(nChars)) {
    col <- charmat[, j]
    prep <- prep_tip_states_single(col)
    tip_states <- prep$tip_states
    k <- prep$k
    
    # compute MLE for observed
    mle <- mle_t(edge, Ntip, tip_states, as.numeric(1.0), k, lower, upper, tol)
    t_hat <- mle$t_hat
    score <- t_hat * nEdge
    
    # compute bestScore: character induced by bestSplit (two states)
    best_ch_int <- ifelse(as.logical(bestSplit), 1L, 0L)
    best_prep <- prep_tip_states_from_int(best_ch_int)
    best_mle <- mle_t(edge, Ntip, best_prep$tip_states, as.numeric(1.0), 2L, lower, upper, tol)
    bestScore <- best_mle$t_hat * nEdge
    
    # resampling: permute observed tokens preserving counts until desired precision
    tokens <- as.character(col)
    counts <- table(tokens)
    # make vector of tokens length Ntip
    tokenVec <- unname(tokens)
    
    rs_vals <- numeric()
    rs_sum <- 0
    rs_sqsum <- 0
    n <- 0
    while (TRUE) {
      perm <- sample(tokenVec, size = length(tokenVec), replace = FALSE)
      prep_r <- prep_tip_states_single(perm)
      mle_r <- mle_t(edge, Ntip, prep_r$tip_states, as.numeric(1.0), prep_r$k, lower, upper, tol)
      rmScore <- mle_r$t_hat * nEdge
      n <- n + 1
      rs_vals[n] <- rmScore
      rs_sum <- rs_sum + rmScore
      rs_sqsum <- rs_sqsum + rmScore * rmScore
      if (n >= 2) {
        mean_rm <- rs_sum / n
        se_rm <- sqrt((rs_sqsum - n * mean_rm * mean_rm) / (n - 1)) / sqrt(n)
        if (se_rm <= precision || n >= maxResample) break
      }
    }
    rmMean <- mean(rs_vals)
    rmSE <- if (length(rs_vals) >= 2) sqrt(var(rs_vals) / length(rs_vals)) else NA_real_
    
    mlci <- (score - rmMean) / (bestScore - rmMean)
    
    res$t_hat[j] <- t_hat
    res$score[j] <- score
    res$rmMean[j] <- rmMean
    res$rmSE[j] <- rmSE
    res$bestScore[j] <- bestScore
    res$mlci[j] <- mlci
  }
  res
}

# Convert a single column (vector of tip tokens) into the tip_states matrix and weights
prep_tip_states_single <- function(col) {
  # col: character or factor or integer vector of length Ntip
  col <- as.character(col)
  Ntip <- length(col)
  # map unique tokens to integers 0..k-1
  lvl <- sort(unique(col))
  k <- length(lvl)
  map <- setNames(seq_len(k) - 1L, lvl)
  mat <- matrix(map[col], nrow = Ntip, ncol = 1)
  weights <- 1.0
  list(tip_states = mat, weights = weights, k = k)
}

# Convert an integer vector (0..k-1) of length Ntip into tip_states matrix (Ntip x 1)
prep_tip_states_from_int <- function(intvec) {
  Ntip <- length(intvec)
  mat <- matrix(as.integer(intvec), nrow = Ntip, ncol = 1)
  list(tip_states = mat, weights = 1.0, k = length(unique(intvec)))
}

# Compute best split character for a tree: iterate over splits and choose the one with minimal 'imbalance'
# We implement the same logic you used: require TreeTools for as.Splits and SplitImbalance
best_split_character <- function(tree) {
  splits <- TreeTools::as.Splits(tree)
  im <- TreeTools::SplitImbalance(splits)
  best <- which.min(im)
  return(splits[[best]]) # logical vector of length Ntip
}
