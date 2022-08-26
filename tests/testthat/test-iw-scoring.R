test_that("IW Scoring", {
  library('TreeTools', quietly = TRUE)
  data('Lobo', package = 'TreeTools')
  dataset <- Lobo.phy
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
  
  unlisted <- unlist(dataset, use.names = FALSE)
  tokenMatrix <- matrix(simpleCont[unlisted], nChar, nTip)
  charInfo <- apply(tokenMatrix, 1, CharacterInformation)
  needsInapp <- rowSums(tokenMatrix == '-') > 2
  inappSlowdown <- 3L # A guess
  rawPriority <- charInfo / ifelse(needsInapp, inappSlowdown, 1)
  priority <- startWeights * rawPriority
  informative <- needsInapp | charInfo > 0
  # Will work from end of sequence to start.
  charSeq <- seq_along(charInfo)[informative][order(priority[informative])] - 1L

  
  weight <- startWeights
  minLength <- MinimumLength(dataset, compress = TRUE)
  
  expect_equal(.IWScore(edge, morphyObjects, weight, minLength, concavity),
               morphy_iw(edge, morphyObjects, weight, minLength, charSeq, 
                         concavity, Inf))
  
  expect_equal(Inf, morphy_iw(edge, morphyObjects, weight, minLength, charSeq,
                              concavity, 0))
  
})
