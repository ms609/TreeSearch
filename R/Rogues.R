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

#' Detect rogue taxa using splits contradicted
#' @references 
#' \insertRef{Aberer2013}{TreeSearch}
#' \insertRef{Wilkinson2017}{TreeSearch}
#' @examples 
#' library("TreeTools", quietly = TRUE)
#' trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
#' trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
#' SplitBuster(trees)
#' @importFrom TreeTools SplitInformation
#' @importFrom TreeDist Entropy
#' @export
TipInformation <- function (trees) {
  t1 <- trees[[1]]
  nTip <- NTip(t1)
  halfTip <- nTip / 2L
  
  rawSplits <- do.call(rbind,
                       lapply(trees, function (tr) as.Splits(tr, t1)))
  hash <- as.character(rawSplits)
  dim(hash) <- dim(rawSplits)
  hash <- apply(hash, 1, paste0, collapse = '')
  uniqueSplits <- structure(class = "Splits", nTip = nTip,
                            tip.label = TipLabels(t1),
                            rawSplits[!duplicated(hash), , drop = FALSE])
  weights <- table(hash)
  # TODO improve efficiency by remaining in "raw" mode
  splits <- as.logical(uniqueSplits)
  splitH <- apply(splits, 1, SplitEntropy)['H1', ]
  
  nSplits <- length(splitH)
  nMinus1 <- nSplits - 1L
  
  splitCombs <- combn(nSplits, 2)
  combWeights <- weights[splitCombs]
  dim(combWeights) <- dim(splitCombs)
  combWeights <- apply(combWeights, 2, prod)
  
  miLosses <- apply(splitCombs, 2, function (i) {
    split1 <- splits[i[1], ]
    split2 <- splits[i[2], ]
    
    A1A2 <- sum(split1 & split2)
    A1B2 <- sum(split1 & !split2)
    B1A2 <- sum(!split1 & split2)
    B1B2 <- sum(!split1 & !split2)
    overlaps <- c(A1A2, A1B2, B1A2, B1B2)
    
    A1 <- A1A2 + A1B2
    A2 <- A1A2 + B1A2
    B1 <- B1A2 + B1B2
    B2 <- A1B2 + B1B2
    
    jointHStart <- Entropy(overlaps / nSplits)
    h1Start <- Entropy(c(A1, B1) / nSplits)
    h2Start <- Entropy(c(A2, B2) / nSplits)
    miStart <- h1Start + h2Start - jointHStart
    
    jointHAA <- Entropy((overlaps - c(1, 0, 0, 0)) / nMinus1)
    jointHAB <- Entropy((overlaps - c(0, 1, 0, 0)) / nMinus1)
    jointHBA <- Entropy((overlaps - c(0, 0, 1, 0)) / nMinus1)
    jointHBB <- Entropy((overlaps - c(0, 0, 0, 1)) / nMinus1)
    
    h1AA <- Entropy(c(A1 - 1L, B1) / nMinus1)
    h1AB <- Entropy(c(A1 - 1L, B1) / nMinus1)
    h1BA <- Entropy(c(A1, B1 - 1L) / nMinus1)
    h1BB <- Entropy(c(A1, B1 - 1L) / nMinus1)
    
    h2AA <- Entropy(c(A2 - 1L, B2) / nMinus1)
    h2AB <- Entropy(c(A2 - 1L, B2) / nMinus1)
    h2BA <- Entropy(c(A2, B2 - 1L) / nMinus1)
    h2BB <- Entropy(c(A2, B2 - 1L) / nMinus1)
    
    miAA <- h1AA + h2AA - jointHAA
    miAB <- h1AB + h2AB - jointHAB
    miBA <- h1BA + h2BA - jointHBA
    miBB <- h1BB + h2BB - jointHBB
    
    viAA <- jointHAA - miAA
    viAB <- jointHAB - miAB
    viBA <- jointHBA - miBA
    viBB <- jointHBB - miBB
    
    miWithout <- ifelse(split1, ifelse(split2, miAA, miAB),
                        ifelse(split2, miBA, miBB))
    miLost <- miStart - miWithout
    
    # Return:
    miLost
  })
  colSums(t(miLosses) * combWeights)
  
}
