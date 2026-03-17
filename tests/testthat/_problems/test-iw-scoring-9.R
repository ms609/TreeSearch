# Extracted from test-iw-scoring.R:9

# test -------------------------------------------------------------------------
library("TreeTools", quietly = TRUE)
data("Lobo", package = "TreeTools")
dataset <- Lobo.phy
tree <- NJTree(dataset)
.IWScore <- function (edge, morphyObjs, weight, minLength, concavity) {
    steps <- preorder_morphy_by_char(edge, morphyObjs)
    homoplasies <- steps - minLength
    fit <- homoplasies / (homoplasies + concavity)
    sum(fit * weight)
  }
