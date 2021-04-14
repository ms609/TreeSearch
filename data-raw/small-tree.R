library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
devtools::load_all("c:/research/r/TreeSearch")
#library("TreeSearch")
library("TreeDist")
library("Quartet", exclude = 'RobinsonFoulds')
cols <- paste0(Ternary::cbPalette8, '44')[rep(1:8,
                                              c(2, 4, 3, rep(1, 5)))]

nTree <- NUnrooted(nTip)
message(Sys.time(), ": Generating all ", nTip, "-leaf trees...")
trees <- lapply(seq_len(nTree) - 1L, as.phylo, nTip = nTip)

CompareMethods <- function (repl, nTip,
                            nTree = NUnrooted(nTip),
                            trees = lapply(seq_len(nTree) - 1L, as.phylo, 
                                           nTip = nTip)) {
  cacheFile <- paste0('gen-', nTip, 'tip/', repl, '.csv')
  if (file.exists(cacheFile)) {
    tab <- read.table(cacheFile, sep = ',')
    performance <- matrix(unlist(tab), nrow = nrow(tab), ncol = ncol(tab),
                          dimnames = dimnames(tab))
  } else {
    message(Sys.time(), ": Simulating data from ", nTip, "-leaf tree ", repl)
    generative <- RandomTree(nTip, TRUE)
    nEdge <- dim(generative$edge)[1] - 1L # Root edge counted twice
    # Desire mean 1 change per character.
    edgeLengths <- rgamma(nEdge, shape = 1) / nEdge
    generative$edge.length <- c(0, edgeLengths)
    dirName <- paste0('gen-', nTip, 'tip')
    if (!dir.exists(dirName)) dir.create(dirName)
    write.tree(generative, paste0(dirName, '/', repl, '.tre'))
    
    
    nChar <- 2 * nEdge
    #shape = 2.5 from Iotuba IQTree ML analysis
    dataBits <- lapply(phangorn::discrete.gamma(2.5, 4), function (rate)
      phangorn::simSeq(generative, l = nChar, rate = rate,
                       type = 'USER', levels = 0:1, rootseq = rep(0, nChar)))
    dataset <- do.call(c, dataBits)
    write.nexus.data(dataset, file = paste0(dirName, '/', repl, '.nex'))
    
    morphyObj <- PhyDat2Morphy(dataset)
    on.exit(morphyObj <- UnloadMorphy(morphyObj))
    
    at <- attributes(dataset)
    characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                              byTaxon = FALSE, concatenate = FALSE)
    weight <- at$weight
    morphyObjects <- lapply(characters, SingleCharMorphy)
    on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)))
      
    nLevel <- length(at$level)
    nChar <- at$nr
    cont <- at$contrast
    simpleCont <- ifelse(rowSums(cont) == 1,
                         apply(cont != 0, 1, function (x) at$levels[x][1]),
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
    
    tokenMatrix <- matrix(simpleCont[unlisted], nChar, nTip, byrow = FALSE)
    profileTables <- apply(tokenMatrix, 1, table)
    data('profiles', package = 'TreeSearch')
    profileCost <- lapply(profileTables, function (x) {
      x <- sort(x[x > 1])
      n <- length(x)
      prof <- switch(n,
        0,
        profiles[[sum(x)]][[n]][[x[1] - 1L]]
      )
      prof <- prof - prof[1]
    })
    
    
    PP <- function (edge) {
      TreeSearch:::morphy_pp(edge, morphyObjects, weight, profileCost, charSeq, Inf)
    }
    
    profileMax <- lapply(profileCost, function (x) {
      if (length(x) == 1L) x else x / max(x)
    })
    PPx <- function (edge) {
      TreeSearch:::morphy_pp(edge, morphyObjects, weight, profileMax, charSeq, Inf)
    }
    
    ewMax <- lapply(profileTables, function (x) {
      switch(length(x),
             0,
             seq_len(min(x)) / min(x))
    })
    EWMax <- function (edge) {
      TreeSearch:::morphy_pp(edge, morphyObjects, weight, ewMax, charSeq, Inf)
    }
    
    
    # Initialize variables and prepare search
    
    message(Sys.time(), ": Scoring each ", nTip, "-leaf tree")
    scores <- vapply(trees, function (tr) {
      edge <- tr$edge
      c(px = PPx(edge),
        pp = PP(edge),
        i1 = IW(edge, 1),
        i3 = IW(edge, 3),
        i10.5 = IW(edge, 10.5),
        i36 = IW(edge, 36),
        ew = TreeSearch:::preorder_morphy(edge, morphyObj),
        ex = EWMax(edge)
      )
    }, c(px = 0, pp = 0, i1 = 0, i3 = 0, i10.5 = 0, i36 = 0, ew = 0, ex = 0))
    
    
    message(Sys.time(), ": Calculating distances: CID")
    cid <- ClusteringInfoDistance(trees, generative, normalize = TRUE)
    # message(Sys.time(), ": Calculating distances: QD")
    # qd <- QuartetDivergence(SingleTreeQuartetAgreement(trees, generative),
    #                         similarity = FALSE)
    qd = replicate(0, nTrees)
    tbr = replicate(0, nTrees)
    # message(Sys.time(), ": Calculating distances: TBR")
    # tbr <- vapply(trees, TBRDist::TBRDist, 0L,
    #               tree2 = generative, exact = TRUE)
    # 
    # Plot <- function (x, lab) {
    #   Normalize <- function (x) {
    #     extra <- x - min(x)
    #     extra / max(extra)
    #   }
    #   score <- apply(scores, 1, Normalize)
    #   plot(0, 0, type = 'n', xlab = lab, ylab = 'Excess score',
    #        xlim = c(0, 0.8), ylim = c(0.3, 0))
    #   points(score[, 'i1'] ~ x, pch = 0, col = cols[1])
    #   points(score[, 'i3'] ~ x, pch = 1, col = cols[2])
    #   points(score[, 'i10.5'] ~ x, pch = 2, col = cols[3])
    #   points(score[, 'i36'] ~ x, pch = 3, col = cols[4])
    #   points(score[, 'ew'] ~ x, pch = 4, col = cols[5])
    #   legend('topright', bty = 'n',
    #          col = cols[1:5], pch = 16, 
    #          legend = c('IW, k = 1', 'IW, k = 3', 'I", k = 10.5', 'IW, k = 36', 'EW'))
    # }
    # 
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
        equalToGen = sum(iScore == generativeScore[i]),
        mean(cid[minima]), 
        mean(qd[minima]),
        mean(tbr[minima])
      )
    }, c('betterThanGen' = 0, 'equalToGen' = 0,
         'cidFromGen' = 0, 'qdFromGen' = 0, 'tbrFromGen' = 0))
    colnames(performance) <- rownames(scores)
    message(Sys.time(), ": Complete.\n")
    write.table(performance, file = cacheFile, sep = ',')
  }
  performance
}
compare.VAL <- matrix(0, 5, 8)
