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
#' 
#' TipInformation(trees)
#' @importFrom TreeDist SplitEntropy Entropy
#' @export
TipInformation <- function (trees) {
  t1 <- trees[[1]]
  nTip <- NTip(t1)
  nMinus1 <- nTip - 1L
  
  rawSplits <- do.call(rbind,
                       lapply(trees, function (tr) as.Splits(tr, t1)))
  hash <- as.character(rawSplits)
  dim(hash) <- dim(rawSplits)
  hash <- apply(hash, 1, paste0, collapse = '')
  uniqueSplits <- structure(class = "Splits", nTip = nTip,
                            tip.label = TipLabels(t1),
                            rawSplits[!duplicated(hash), , drop = FALSE])
  weights <- table(hash, dnn = NULL)[unique(hash)]
  # TODO improve efficiency by remaining in "raw" mode
  splits <- as.logical(uniqueSplits)
  splitH <- apply(splits, 1, SplitEntropy)['H1', ]
  
  nSplits <- length(splitH)
  
  splitCombs <- combn(nSplits, 2)
  combWeights <- weights[splitCombs]
  dim(combWeights) <- dim(splitCombs)
  combWeights <- apply(combWeights, 2, prod)
  
  fromDiff <- apply(splitCombs, 2, function (i) {
    split1 <- splits[i[1], ]
    split2 <- splits[i[2], ]
    rbind(split1, split2)
    ifelse(split1, ifelse(split2, "TT", "TF"), ifelse(split2, "FT", "FF"))
    
    T1T2 <- sum(split1 & split2)
    T1F2 <- sum(split1 & !split2)
    F1T2 <- sum(!split1 & split2)
    F1F2 <- sum(!split1 & !split2)
    overlaps <- c(T1T2, T1F2, F1T2, F1F2)
    names(overlaps) <- c('TT', 'TF', 'FT', 'FF') #TODO delete
    
    q <- overlaps / nTip
    p <- cbind(overlaps - c(1, 0, 0, 0),
               overlaps - c(0, 1, 0, 0),
               overlaps - c(0, 0, 1, 0),
               overlaps - c(0, 0, 0, 1)) / nMinus1
    q[q == 0] <- NA
    log2pq <- log2(p / q)
    log2pq[!is.finite(log2pq)] <- 0
    dKL <- p * log2pq
    
    dKLPQ <- dKL[ifelse(split1, ifelse(split2, 1, 2),
                          ifelse(split2, 3, 4)), ]
    
    rowSums(dKLPQ)
  })
  
  sameSplits <- weights > 1L
  fromSame <- apply(splits[sameSplits, , drop = FALSE], 1, function (split) {
    m <- sum(split)
    n <- nTip - m
    q <- c(m, n) / nTip
    p <- cbind(c(m - 1L, n), c(m, n - 1L)) / nMinus1
    
    dKL <- p* log2(p / q)
    
    dKLPQ <- dKL[ifelse(split, 1, 2), ]
    
    # Return:
    rowSums(dKLPQ)
  })
  
  colSums(t(fromSame) * as.integer(weights[sameSplits] - 1L)) +
    colSums(t(fromDiff) * combWeights)
  
}
