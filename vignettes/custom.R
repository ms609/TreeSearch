## ----init, message = FALSE----------------------------------------------------
library("TreeTools", quietly = TRUE, warn.conflict = FALSE)
library("TreeSearch")

# Plot trees nicely
PlotTree <- function(tree, ...) {
  oPar <- par(mar = rep(0, 4), cex = 0.9)
  plot(tree)
  par(oPar)
}

## ----tci-setup----------------------------------------------------------------
tree <- PectinateTree(8)
PlotTree(tree)

TCIScore <- function(parent, child, dataset) {
  tree$edge <- cbind(parent, child)
  TotalCopheneticIndex(tree)
}

TCIScore(tree$edge[, 1], tree$edge[, 2], NA)

## ----tci-search---------------------------------------------------------------
result <- TreeSearch(tree, dataset = EmptyPhyDat(tree),
                     InitializeData = DoNothing, CleanUpData = DoNothing,
                     TreeScorer = TCIScore,
                     maxIter = 50L, maxHits = 10L, 
                     verbosity = 1L)

PlotTree(result)

## ----cid----------------------------------------------------------------------
startTree <- BalancedTree(8)

DistanceScore <- function(parent, child, dataset) {
  tmpTree <- startTree
  tmpTree$edge <- cbind(parent, child)
  distance <- TreeDist::ClusteringInfoDistance(startTree, tmpTree)
  # Return:
  -distance
}

result <- TreeSearch(RandomTree(8, root = TRUE), dataset = EmptyPhyDat(tree),
                     InitializeData = DoNothing, CleanUpData = DoNothing,
                     TreeScorer = DistanceScore,
                     maxIter = 50L, maxHits = 10L, 
                     verbosity = 1L)

par(mfrow = c(1, 2))
PlotTree(startTree)
PlotTree(result)


## ----iw-setup-----------------------------------------------------------------
IWInitMorphy <- function (dataset) {
  attr(dataset, "morphyObjs") <- 
    lapply(PhyToString(dataset, byTaxon = FALSE, useIndex = FALSE, 
                       concatenate = FALSE), 
           SingleCharMorphy)
  
  # Return:
  dataset
}

## ----iw-destroy---------------------------------------------------------------
IWDestroyMorphy <- function (dataset) {
  vapply(attr(dataset, "morphyObjs"), UnloadMorphy, integer(1))
}

## ----iw-score-----------------------------------------------------------------
IWScoreMorphy <- function (parent, child, dataset, concavity = 10L, 
                           minLength = attr(dataset, "min.length"), ...) {
  steps <- vapply(attr(dataset, "morphyObjs"), MorphyLength,
                  parent = parent, child = child, integer(1))
  homoplasies <- steps - minLength
  fit <- homoplasies / (homoplasies + concavity)
  # Return:
  sum(fit * attr(dataset, "weight"))
}

## ----iw-search, message = FALSE-----------------------------------------------
data("inapplicable.datasets")
dataset <- congreveLamsdellMatrices[[42]]

# Populate `min.length` attribute
dataset <- PrepareDataIW(dataset)
iwTree <- TreeSearch(NJTree(dataset), dataset,
                     InitializeData = IWInitMorphy,
                     CleanUpData = IWDestroyMorphy,
                     TreeScorer = IWScoreMorphy,
                     concavity = 10, # Will be sent to TreeScorer
                     verbosity = 1)


## ----iw-bootstrap-------------------------------------------------------------
IWBootstrap <- function (edgeList, dataset, concavity = 10L, EdgeSwapper = NNISwap, 
                         maxIter, maxHits, verbosity = 1L, ...) {
  att <- attributes(dataset)
  startWeights <- att[["weight"]]
  
  # Decompress phyDat object so each character is listed once
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep.int(eachChar, startWeights)
  
  # Resample characters
  resampling <- tabulate(sample(deindexedChars, replace = TRUE), length(startWeights))
  sampled <- resampling != 0
  sampledData <- lapply(dataset, function (x) x[sampled])
  sampledAtt <- att
  sampledAtt[["index"]] <- rep.int(seq_len(sum(sampled)), resampling[sampled])
  sampledAtt[["weight"]] <- resampling[sampled]
  sampledAtt[["nr"]] <- length(sampledAtt[["weight"]])
  sampledAtt[["min.length"]] <- minLength <- att[["min.length"]][sampled]
  sampledAtt[["morphyObjs"]] <- att[["morphyObjs"]][sampled]
  attributes(sampledData) <- sampledAtt
  
  # Search using resampled dataset
  res <- EdgeListSearch(edgeList[1:2], sampledData, TreeScorer = IWScoreMorphy,
                        concavity = concavity, minLength = minLength,
                        EdgeSwapper = EdgeSwapper, 
                        maxIter = maxIter, maxHits = maxHits,
                        verbosity = verbosity - 1L)
  
  res[1:2]
}


## ----iw-ratchet, message = FALSE----------------------------------------------
ratchetTree <- Ratchet(tree = iwTree, dataset = dataset,
                       concavity = 10,
                       InitializeData = IWInitMorphy, 
                       CleanUpData = IWDestroyMorphy,
                       TreeScorer = IWScoreMorphy,
                       Bootstrapper = IWBootstrap,
                       ratchIter = 2, ratchHits = 2,
                       searchIter = 20, searchHits = 10,
                       verbosity = 2)


