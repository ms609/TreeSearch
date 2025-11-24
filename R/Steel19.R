#' Steel-inspired consistency index
#' @examples
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' tree <- TreeTools::NJTree(dataset)
#' SteelInfo(tree, dataset)
#' @export
SteelInfo <- function(tree, dataset) {
  lengths <- CharacterLength(tree, dataset)
  minLengths <- MinimumLength(dataset)
  # nTokens <- lengths(apply(PhyDatToMatrix(dataset), 2, intersect, 0:9))
  
  lengthTab <- table(lengths, minLengths)

  lapply(colnames(lengthTab), function(steps) {
    nTokens <- as.numeric(steps) + 1
    tokenP <- rep(1 / nTokens, nTokens)
    col <- lengthTab[, steps]
    maxSteps <- NTip(tree) - nTokens + 1
    lnP <- active_parsimony_dist(
      tree, tokenP,
      maxSteps)
      #as.numeric(rownames(lengthTab))[max(which(col > 0))])
    expect <- parsimony_moments(tree, tokenP)
    normalized <- (as.numeric(names(lnP)) - expect$expectation) / sqrt(expect$variance)
  })
  # active_parsimony_dist(tree, 
}
