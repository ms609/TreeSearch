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
  dists <- .TipDistances(trees)
  
  means <- rowMeans(dists, dims = 2)
  devs <- apply(dists, 1:2, function(x) mad(x))
  diag(devs) <- NA
  relDevs <- devs / mean(means[lower.tri(means)])
  rowMeans(relDevs, na.rm = TRUE)
}

.TipDistances <- function (trees) {
  nTip <- NTip(trees)
  if (length(unique(nTip)) > 1) {
    stop("Trees must have same number of leaves")
  }
  nTip <- nTip[1]
  trees[-1] <- lapply(trees[-1], RenumberTips, trees[[1]])
  dists <- vapply(trees, Cophenetic, matrix(0, nTip, nTip))
}

#' @importFrom grDevices hcl
.TipCols <- function (trees, luminence = 50) {
  dists <- .TipDistances(trees)
  
  means <- rowMeans(dists, dims = 2)
  devs <- apply(dists, 1:2, function(x) mad(x))
  diag(devs) <- NA
  relDevs <- devs / mean(means[lower.tri(means)])
  
  pc <- cmdscale(means, k = 1)
  pc <- pc - min(pc)
  pc <- pc * 340 / max(pc)
  
  setNames(hcl(h =  pc, c = 100 * (1 - rowMeans(relDevs, na.rm = TRUE)),
               l = luminence),
           TipLabels(trees[[1]]))
  
}

#' Detect rogue taxa using phylogenetic information distance
#' @references 
#' \insertRef{Aberer2013}{TreeSearch}
#' \insertRef{Wilkinson2017}{TreeSearch}
#' @examples 
#' library("TreeTools", quietly = TRUE)
#' trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
#' trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
#' 
#' TipInformation(trees)
#' sb <- TipInformation(trees)
#' col <- hcl.colors(ceiling(max(sb) *1.13), 'inferno')[ceiling(sb)]
#' plot(consensus(trees), tip.color = col)
#' plot(ConsensusWithout(trees, names(sb[sb == max(sb)])))
#' @importFrom TreeDist PhylogeneticInfoDistance CladisticInfo
#' @export
TipRemover <- function (trees) {
  tips <- trees[[1]]$tip.label
  startInfo <- mean(CladisticInfo(trees))
  info <- vapply(tips, function (drop) {
    tr <- lapply(trees, drop.tip, drop)
    c(meanInfo = mean(CladisticInfo(tr)),
      meanDist = mean(PhylogeneticInfoDistance(tr)))
  }, double(2))
  (mean(PhylogeneticInfoDistance(trees)) - info['meanDist', ]) /
  (startInfo / info['meanInfo', ])
}

#' Calculate the most informative consensus tree
#' 
#' Uses the splitwise information content as a shortcut, which involves double
#' counting of some information (which may or may not be desirable) but is
#' calculable in polynomial ratehr than exponential time.
#' 
#' @template MRS
#' @importFrom ape consensus
#' @importFrom TreeDist SplitwiseInfo SplitFrequency
#' @export
BestConsensus <- function (trees) {
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return (trees)
    }
    if (!is.list(trees)) {
      stop("`trees` must be a list of `phylo` objects")
    }
    trees <- structure(trees, class = 'multiPhylo')
  }
  trees <- lapply(trees, RenumberTips, trees[[1]])
  trees <- lapply(trees, Preorder)
  nTree <- length(trees)
  
  cons <- vector('list', nTip - 2L)
  info <- double(nTip - 2)
  cons[[1]] <- consensus(trees, p = 0.50)
  info[1] <- SplitwiseInfo(cons[[1]], SplitFrequency(cons[[1]], trees) / nTree)
  for (i in 1 + seq_len(nTip - 3L)) {
    tipScores <- TipRemover(trees)
    candidate <- which.max(tipScores)
    trees <- lapply(trees, drop.tip, candidate)
    cons[[i]] <- tr <- consensus(trees, p  = 0.50)
    info[i] <- SplitwiseInfo(tr, SplitFrequency(tr, trees) / nTree)
  }
  cons[[which.max(info)]]
}
  