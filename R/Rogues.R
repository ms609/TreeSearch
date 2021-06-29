#' @importFrom stats cophenetic
Cophenetic <- function (x) {
  if (is.null(x$edge.length)) {
    x$edge.length <- rep_len(1, dim(x$edge)[1])
  }
  ret <- cophenetic(x)
  ret[ret < sqrt(.Machine$double.eps)] <- 0
  ret
}

# @examples 
# library("TreeTools", quietly = TRUE)
# trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
# plot(consensus(trees))
# instab <- TipInstability(trees)
# plot(ConsensusWithout(trees, names(instab[instab > 0.2])))
# @template MRS
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
#' 
#' Calculate the volatility of each tip: namely, the impact on the mean
#' phylogenetic information distance between trees when that tip is removed.
#' 
#' @return `TipVolatility()` returns a named vector listing the volatility
#' index calculated for each leaf.
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
#' @importFrom TreeDist PhylogeneticInfoDistance 
#' @importFrom TreeTools CladisticInfo
#' @export
TipVolatility <- function (trees) {
  tips <- trees[[1]]$tip.label
  startInfo <- mean(CladisticInfo(trees))
  info <- vapply(tips, function (drop) {
    tr <- lapply(trees, DropTip, drop)
    c(meanInfo = mean(CladisticInfo(tr)),
      meanDist = mean(PhylogeneticInfoDistance(tr, normalize = TRUE)))
  }, double(2))
  mean(PhylogeneticInfoDistance(trees, normalize = TRUE)) - info['meanDist', ]
}

#' Calculate the most informative consensus tree
#' 
#' Uses the splitwise information content as a shortcut, which involves double
#' counting of some information (which may or may not be desirable) but is
#' calculable in polynomial ratehr than exponential time.
#' 
#' @template MRS
#' @importFrom ape consensus
#' @importFrom TreeDist ConsensusInfo 
#' @importFrom TreeTools SplitFrequency
#' @export
BestConsensus <- function (trees, info = 'clustering') {
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return (trees)
    }
    if (!is.list(trees)) {
      stop("`trees` must be a list of `phylo` objects")
    }
    trees <- structure(trees, class = 'multiPhylo')
  }
  lastMessage <- Sys.time()
  trees <- lapply(trees, RenumberTips, trees[[1]])
  trees <- lapply(trees, Preorder)
  nTip <- NTip(trees[[1]])
  nTree <- length(trees)
  
  tr <- trees
  candidates <- character(nTip - 2L)
  score <- double(nTip - 2)
  score[1] <- ConsensusInfo(trees, info = info)
  for (i in 1 + seq_len(nTip - 3L)) {
    if (difftime(Sys.time(), lastMessage) > 2) {
      lastMessage <- Sys.time()
      message(lastMessage, ": Removing tip ", i)
    }
    tipScores <- TipVolatility(tr)
    candidate <- which.max(tipScores)
    if (length(candidate)) {
      candidates[i] <- names(candidate)
    }
    tr <- lapply(tr, DropTip, candidate)
    score[i] <- ConsensusInfo(tr, info = info)
  }
  droppers <- candidates[seq_len(which.max(score))[-1]]
  consensus(lapply(trees, DropTip, droppers), p = 0.5)
}

#' Calculate the most informative consensus tree
#' 
#' Uses the splitwise information content as a shortcut, which involves double
#' counting of some information (which may or may not be desirable) but is
#' calculable in polynomial rather than exponential time.
#' 
#' @examples 
#' library("TreeTools", warn.conflicts = FALSE)
#' trees <- list(
#'      read.tree(text = '((a, y), (b, (c, (z, ((d, e), (f, (g, x)))))));'),
#'      read.tree(text = '(a, (b, (c, (z, (((d, y), e), (f, (g, x)))))));'),
#'      read.tree(text = '(a, (b, ((c, z), ((d, (e, y)), ((f, x), g)))));'),
#'      read.tree(text = '(a, (b, ((c, z), ((d, (e, x)), (f, (g, y))))));'),
#'      read.tree(text = '(a, ((b, x), ((c, z), ((d, e), (f, (g, y))))));')
#'      )
#' cons <- consensus(trees, p = 0.5)
#' plot(cons)
#' LabelSplits(cons, SplitFrequency(cons, trees) / length(trees))
#' reduced <- Roguehalla(trees, info = 'phylogenetic')
#' plot(reduced)
#' LabelSplits(reduced, SplitFrequency(reduced, trees) / length(trees))
#' @template MRS
#' @importFrom TreeDist ConsensusInfo
#' @importFrom TreeTools DropTip SplitFrequency Preorder RenumberTips
#' @export
Roguehalla <- function (trees, dropset = 1, info = 'clustering') {
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
  majority <- 0.5 + sqrt(.Machine$double.eps)
  
  startTip <- NTip(trees[[1]])
  best <- ConsensusInfo(trees, info = info)
  # cons <- consensus(trees, p = 0.5)
  # stopifnot(SplitwiseInfo(cons, SplitFrequency(cons, trees) / nTree) == best)
  
  .Drop <- function (n) {
    drops <- combn(NTip(trees[[1]]), n)
    candidates <- apply(drops, 2, function (drop) {
      ConsensusInfo(lapply(trees, DropTip, drop), info = info)
    })
    if (max(candidates) > best) {
      list(info = max(candidates), drop = drops[, which.max(candidates)])
    } else {
      NULL
    }
  }
  
  lastMessage <- Sys.time()
  .Message <- function (...) {
    if (difftime(Sys.time(), lastMessage) > 2) {
      message(paste0(Sys.time(), ': ', ...))
      lastMessage <- Sys.time()
    }
    lastMessage
  }
  repeat {
    improved <- FALSE
    for (i in seq_len(dropset)) {
      lastMessage <- .Message("Dropped ", startTip - NTip(trees[[1]]),
                              " leaves to render ", signif(best),
                              " bits; dropping ", i, " more.")
      dropped <- .Drop(i)
      if (!is.null(dropped)) {
        improved <- TRUE
        best <- dropped$info
        trees <- lapply(trees, DropTip, dropped$drop)
        break
      }
    }
    if (!improved) break
  }
  
  # Return:
  consensus(trees, p = majority)
}
