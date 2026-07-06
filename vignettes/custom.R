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
                     TreeScorer = DistanceScore,
                     maxIter = 50L, maxHits = 10L,
                     verbosity = 1L)

par(mfrow = c(1, 2))
PlotTree(startTree)
PlotTree(result)



## ----iw-score, message = FALSE------------------------------------------------
data("inapplicable.datasets")
dataset <- congreveLamsdellMatrices[[42]]

IWScore <- function (parent, child, dataset, concavity = 10,
                     minLength = MinimumLength(dataset, compress = TRUE)) {
  tree <- structure(list(edge = cbind(parent, child),
                        tip.label = names(dataset),
                        Nnode = length(dataset) - 1L), class = "phylo")
  homoplasy <- CharacterLength(tree, dataset, compress = TRUE) - minLength
  fit <- homoplasy / (homoplasy + concavity)
  # Return:
  sum(fit * attr(dataset, "weight"))
}

iwTree <- TreeSearch(NJTree(dataset), dataset, TreeScorer = IWScore,
                     maxIter = 50L, maxHits = 10L, verbosity = 1)


## ----iw-bootstrap-------------------------------------------------------------
IWBootstrap <- function (edgeList, dataset, EdgeSwapper = NNISwap,
                         maxIter, maxHits, verbosity = 1L, ...) {
  startWeights <- attr(dataset, "weight")
  resampling <- tabulate(sample(rep.int(seq_along(startWeights), startWeights),
                                replace = TRUE), length(startWeights))
  # R copy-on-modify: the caller's `dataset` is unchanged.
  attr(dataset, "weight") <- as.integer(resampling)

  res <- EdgeListSearch(edgeList[1:2], dataset, TreeScorer = IWScore,
                        EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        verbosity = verbosity - 1L, ...)
  # Return:
  res[1:2]
}


## ----iw-ratchet, message = FALSE----------------------------------------------
ratchetTree <- Ratchet(iwTree, dataset, TreeScorer = IWScore,
                       Bootstrapper = IWBootstrap,
                       ratchIter = 2, ratchHits = 2,
                       searchIter = 20, searchHits = 10,
                       verbosity = 2)

