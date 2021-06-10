#' @importFrom ape cophenetic
Cophenetic <- function (x) {
  if (is.null(x$edge.length)) {
    x$edge.length <- rep_len(1, dim(x$edge)[1])
  }
  ret <- cophenetic(x)
  ret[ret < sqrt(.Machine$double.eps)] <- 0
  ret
}

#' @examples 
#' library("TreeTools", quietly = TRUE)
#' trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
#' plot(consensus(trees))
#' instab <- TipInstability(trees)
#' plot(ConsensusWithout(trees, names(instab[instab > 0.2])))
#' @template MRS
TipInstability <- function (trees) {
  nTip <- NTip(trees)
  if (length(unique(nTip)) > 1) {
    stop("Trees must have same number of leaves")
  }
  nTip <- nTip[1]
  trees[-1] <- lapply(trees[-1], RenumberTips, trees[[1]])
  dists <- vapply(trees, Cophenetic, matrix(0, nTip, nTip))
  
  means <- rowMeans(dists, dims = 2)
  devs <- apply(dists, 1:2, function(x) mad(x))
  diag(devs) <- NA
  relDevs <- devs / mean(means[lower.tri(means)])
  rowMeans(relDevs, na.rm = TRUE)
}
