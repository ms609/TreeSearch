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
  
  viLoss <- apply(splitCombs, 2, function (i) {
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
    
    T1 <- T1T2 + T1F2
    T2 <- T1T2 + F1T2
    F1 <- F1T2 + F1F2
    F2 <- T1F2 + F1F2
    
    jointHStart <- Entropy(overlaps / nTip)
    h1Start <- Entropy(c(T1, F1) / nTip)
    h2Start <- Entropy(c(T2, F2) / nTip)
    miStart <- h1Start + h2Start - jointHStart
    viStart <- jointHStart - miStart
    
    
    jointHTT <- Entropy((overlaps - c(1, 0, 0, 0)) / nMinus1)
    jointHTF <- Entropy((overlaps - c(0, 1, 0, 0)) / nMinus1)
    jointHFT <- Entropy((overlaps - c(0, 0, 1, 0)) / nMinus1)
    jointHFF <- Entropy((overlaps - c(0, 0, 0, 1)) / nMinus1)
    ifelse(split1, ifelse(split2, jointHTT, jointHTF),
           ifelse(split2, jointHFT, jointHFF))
    
    h1TT <- Entropy(c(T1 - 1L, F1) / nMinus1)
    h1TF <- h1TT
    h1FT <- Entropy(c(T1, F1 - 1L) / nMinus1)
    h1FF <- h1FT
    h2TT <- Entropy(c(T2 - 1L, F2) / nMinus1)
    h2FT <- h2TT
    h2TF <- Entropy(c(T2, F2 - 1L) / nMinus1)
    h2FF <- h2TF
    rbind(h1 = ifelse(split1, ifelse(split2, h1TT, h1TF), ifelse(split2, h1FT, h1FF)),
          h2 = ifelse(split1, ifelse(split2, h2TT, h2TF), ifelse(split2, h2FT, h2FF)))
    
    
    
    miTT <- h1TT + h2TT - jointHTT
    miTF <- h1TF + h2TF - jointHTF
    miFT <- h1FT + h2FT - jointHFT
    miFF <- h1FF + h2FF - jointHFF
    
    viTT <- jointHTT - miTT
    viTF <- jointHTF - miTF
    viFT <- jointHFT - miFT
    viFF <- jointHFF - miFF
    
    miWithout <- ifelse(split1, ifelse(split2, miTT, miTF),
                        ifelse(split2, miFT, miFF))
    viWithout <- ifelse(split1, ifelse(split2, viTT, viTF),
                        ifelse(split2, viFT, viFF))
    miGained <- miWithout - miStart
    viLost <- viStart - viWithout
    
    # Return:
    miGained + viLost
  })
  
  sameSplits <- weights > 1L
  fromSame <- apply(splits[sameSplits, , drop = FALSE], 1, function (split) {
    m <- sum(split)
    n <- nTip - m
    hStart <- Entropy(c(m, n) / nTip)
    hT <- Entropy(c(m - 1L, n) / nMinus1)
    hF <- Entropy(c(m, n - 1L) / nMinus1)
    
    hWithout <- ifelse(split, hT, hF)
    miGained <- hStart - hWithout
    
    # Return:
    miGained
  })
  
  colSums(t(fromSame) * as.integer(weights[sameSplits] - 1L)) +
    colSums(t(viLoss) * combWeights)
  
  # TODO This overwrites... 
  siGain <- apply(splitCombs, 2, function (i) {
    split1 <- splits[i[1], ]
    split2 <- splits[i[2], ]
    ifelse(split1, ifelse(split2, "TT", "TF"), ifelse(split2, "FT", "FF"))
    
    T1T2 <- sum(split1 & split2)
    T1F2 <- sum(split1 & !split2)
    F1T2 <- sum(!split1 & split2)
    F1F2 <- sum(!split1 & !split2)
    overlaps <- c(T1T2, T1F2, F1T2, F1F2)
    names(overlaps) <- c('TT', 'TF', 'FT', 'FF') #TODO delete
    T1 <- T1T2 + T1F2
    F1 <- F1T2 + F1F2
    T2 <- T1T2 + F1T2
    F2 <- T1F2 + F1F2
    i0 <- SplitInformation(T1, F1) + SplitInformation(T2, F2)
    noOverlap <- !as.logical(overlaps)
    
    iWithout <- 
    if (any(noOverlap)) {
      stopifnot(sum(noOverlap) == 1L) # as splits are non-identical
      overT1 <- c(F2, T2, T2, F2)[noOverlap]
      si0 <- SplitSharedInformation(nTip, T1, overT1)
      
      overT1TT <- c(F2, T2 - 1L, T2 - 1L, F2)[noOverlap]
      overT1TF <- c(F2 - 1L, T2, T2, F2 - 1L)[noOverlap]
      overT1FT <- c(F2, T2 - 1L, T2 - 1L, F2)[noOverlap]
      overT1FF <- c(F2 - 1L, T2, T2, F2 - 1L)[noOverlap]
      
      siTT <- SplitSharedInformation(nMinus1, T1 - 1L, overT1TT)
      siTF <- SplitSharedInformation(nMinus1, T1 - 1L, overT1TF)
      siFT <- SplitSharedInformation(nMinus1, T1, overT1FT)
      siFF <- SplitSharedInformation(nMinus1, T1, overT1FF)
      
      siWithout <- ifelse(split1, ifelse(split2, siTT, siTF),
                          ifelse(split2, siFT, siFF))
    } else {
      si0 <- 0
      si <- double(4)
      for (j in which(overlaps == 1)) {
        si[j] <- switch(j,
                        SplitSharedInformation(nMinus1, T1 - 1L, F2),
                        SplitSharedInformation(nMinus1, T1 - 1L, T2),
                        SplitSharedInformation(nMinus1, T1, T2 - 1L),
                        SplitSharedInformation(nMinus1, T1, F2 - 1L))
      }
      si
      siWithout <- si[ifelse(split1, ifelse(split2, 1, 2),
                             ifelse(split2, 3, 4))]
      
    }
    siGain <- siWithout - si0
    di0 <- i0 - si0 - si0
    di <- c(SplitInformation(T1 - 1L, F1) + SplitInformation(T2 - 1L, F2),
            SplitInformation(T1 - 1L, F1) + SplitInformation(T2, F2 - 1L),
            SplitInformation(T1, F1 - 1L) + SplitInformation(T2 - 1L, F2),
            SplitInformation(T1, F1 - 1L) + SplitInformation(T2, F2 - 1L))
    diWithout <- di[ifelse(split1, ifelse(split2, 1, 2),
                           ifelse(split2, 3, 4))]
    diGain <- diWithout - di0
    
    # Return:
    siGain - diGain
  })
  
  sameSplits <- weights > 1L
  fromSame <- apply(splits[sameSplits, , drop = FALSE], 1, function (split) {
    m <- sum(split)
    n <- nTip - m
    si0 <- SplitInformation(m, n)
    siT <- SplitInformation(m - 1L, n)
    siF <- SplitInformation(m, n - 1L)
    
    siWithout <- ifelse(split, siT, siF)
    i0 <- si0
    iWithout <- siWithout
    di0 <- 0
    diWithout <- 0
    
    siGain <- siWithout - si0
    diGain <- diWithout - di0
    
    # Return:
    siGain - diGain
  })
  colSums(t(fromSame) * as.integer(weights[sameSplits] - 1L)) +
  colSums(t(siGain) * combWeights)
  
}
