library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
library("TreeSearch")
library("TreeDist")
library("Quartet", exclude = 'RobinsonFoulds')
nTip <- 8
cols <- Ternary::cbPalette8

CompareMethods <- function (nTip) {
  message(Sys.time(), ": Simulating data from ", nTip, "-tip tree")
  generative <- RandomTree(nTip, TRUE)
  generative$edge.length <- rgamma(dim(generative$edge)[1], shape = 1)*3
  
  dataBits <- lapply(phangorn::discrete.gamma(1, 4), function (rate)
    phangorn::simSeq(generative, l = 100, rate = rate))
  dataset <- c(dataBits[[1]], dataBits[[2]], dataBits[[3]], dataBits[[4]]) # I can't remember how to do this right
  
  
  morphyObj <- PhyDat2Morphy(dataset)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  
  nTree <- NUnrooted(nTip)
  message(Sys.time(), ": Generating all ", nTip, "-tip trees")
  trees <- lapply(seq_len(nTree) - 1L, as.phylo, nTip = nTip)
  
  
  at <- attributes(dataset)
  characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                            byTaxon = FALSE, concatenate = FALSE)
  weight <- at$weight
  morphyObjects <- lapply(characters, SingleCharMorphy)
  on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)))
    
  nLevel <- length(at$level)
  nChar <- at$nr
  nTip <- length(dataset)
  cont <- at$contrast
  simpleCont <- ifelse(rowSums(cont) == 1,
                       apply(cont != 0, 1, function (x) colnames(cont)[x][1]),
                       '?')
  inappLevel <- at$levels == '-'
  
  if (any(inappLevel)) {
    # TODO this is a workaround until MinimumLength can handle {-, 1}
    cont[cont[, inappLevel] > 0, ] <- 0
    ambiguousToken <- at$allLevels == '?'
    cont[ambiguousToken, ] <- colSums(cont[!ambiguousToken, ]) > 0
  }
  
  powersOf2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
  tmp <- as.integer(cont %*% powersOf2)
  unlisted <- unlist(dataset, use.names = FALSE)
  binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow = FALSE)
  minLength <- apply(binaryMatrix, 1, MinimumLength)
  charSeq <- seq_len(nChar) - 1L
  
  IW <- function (edge, concavity) {
    TreeSearch:::morphy_iw(edge, morphyObjects, weight, minLength, charSeq, concavity, Inf)
  }
  
  # Initialize variables and prepare search
  
  message(Sys.time(), ": Scoring each tree")
  scores <- vapply(trees, function (tr) {
    edge <- tr$edge
    c(ew = TreeSearch:::preorder_morphy(edge, morphyObj),
      i1 = IW(edge, 1),
      i3 = IW(edge, 3),
      i10.5 = IW(edge, 10.5),
      i36 = IW(edge, 36)
    )
  }, c(ew = 0, i1 = 0, i3 = 0, i10.5 = 0, i36 = 0))
  
  
  message(Sys.time(), ": Calculating distances: CID")
  cid <- ClusteringInfoDistance(trees, generative, normalize = TRUE)
  message(Sys.time(), ": Calculating distances: QD")
  qd <- QuartetDivergence(SingleTreeQuartetAgreement(trees, generative),
                          similarity = FALSE)
  message(Sys.time(), ": Calculating distances: TBR")
  tbr <- vapply(trees, TBRDist::TBRDist, 0L,
                tree2 = generative, exact = TRUE)
  
  Plot <- function (x, lab) {
    Normalize <- function (x) {
      extra <- x - min(x)
      extra / max(extra)
    }
    score <- apply(scores, 1, Normalize)
    plot(0, 0, type = 'n', xlab = lab, ylab = 'Excess score',
         xlim = c(0, 0.8), ylim = c(0.3, 0))
    points(score[, 'i1'] ~ x, pch = 0, col = cols[1])
    points(score[, 'i3'] ~ x, pch = 1, col = cols[2])
    points(score[, 'i10.5'] ~ x, pch = 2, col = cols[3])
    points(score[, 'i36'] ~ x, pch = 3, col = cols[4])
    points(score[, 'ew'] ~ x, pch = 4, col = cols[5])
    legend('topright', bty = 'n',
           col = cols[1:5], pch = 16, 
           legend = c('IW, k = 1', 'IW, k = 3', 'I", k = 10.5', 'IW, k = 36', 'EW'))
  }
  
  #par(mfrow = c(2, 1), mar = rep(2, 4), mgp = c(1, 1, 1))
  #Plot(cid, 'CID')
  #Plot(qd, 'QD')
  #Plot(tbr / max(tbr), 'TBR')
  
  message(Sys.time(), ": Evaluating performace")
  generativeScore <- scores[, as.integer(as.TreeNumber(generative)) + 1L]
  performance <- vapply(seq_len(nrow(scores)), function (i) {
    iScore <- scores[i, ]
    minima <- iScore == min(iScore)
    c(betterThanGen = sum(iScore < generativeScore[i]),
      mean(cid[minima]), 
      mean(qd[minima]),
      mean(tbr[minima])
    )
  }, c('betterThanGen' = 0, 'cidFromBest' = 0, 'qdFromBest' = 0, 'tbrFromBest' = 0))
  colnames(performance) <- rownames(scores)
  performance
}

sapply(rep(8, 2), CompareMethods)
message(Sys.time(), ": RUN COMPLETE")
