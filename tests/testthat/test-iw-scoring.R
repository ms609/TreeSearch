test_that("IW Scoring", {
  library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
  data('Lobo', package = 'TreeTools')
  dataset <- Lobo.phy
  
  #dataset <- ReadAsPhyDat('c:/research/r/hyoliths/mbank_X24932_6-19-2018_744.nex')
  tree <- NJTree(dataset)
  
  
  .IWScore <- function (edge, morphyObjs, weight, minLength, concavity) {
    steps <- preorder_morphy_by_char(edge, morphyObjs)
    homoplasies <- steps - minLength
    fit <- homoplasies / (homoplasies + concavity)
    sum(fit * weight)
  }
  
  concavity <- 4.5
  epsilon <- sqrt(.Machine$double.eps)
  
  
  tree <- Preorder(RenumberTips(tree, names(dataset)))
  nTip <- NTip(tree)
  edge <- tree$edge
  
  at <- attributes(dataset)
  characters <- PhyToString(dataset, ps = '', useIndex = FALSE,
                            byTaxon = FALSE, concatenate = FALSE)
  startWeights <- at$weight
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
  
  # Perhaps replace with previous code:
  # inappLevel <- which(at$levels == "-")
  # cont[, inappLevel] <- 0
  
  powersOf2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
  tmp <- as.integer(cont %*% powersOf2)
  unlisted <- unlist(dataset, use.names = FALSE)
  binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow = FALSE)
  minLength <- apply(binaryMatrix, 1, MinimumLength)
  
  tokenMatrix <- matrix(simpleCont[unlisted], nChar, nTip, byrow = FALSE)
  charInfo <- apply(tokenMatrix, 1, CharacterInformation)
  needsInapp <- rowSums(tokenMatrix == '-') > 2
  inappSlowdown <- 3L # A guess
  rawPriority <- charInfo / ifelse(needsInapp, inappSlowdown, 1)
  priority <- startWeights * rawPriority
  informative <- needsInapp | charInfo > 0
  # Will work from end of sequence to start.
  charSeq <- seq_along(charInfo)[informative][order(priority[informative])] - 1L

  
  weight <- startWeights
  
  expect_equal(.IWScore(edge, morphyObjects, weight, minLength, concavity),
               morphy_iw(edge, morphyObjects, weight, minLength, charSeq, 
                         concavity, Inf))
  
  expect_equal(Inf, morphy_iw(edge, morphyObjects, weight, minLength, charSeq,
                              concavity, 0))
  
})