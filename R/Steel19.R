#' Steel-inspired consistency index
#' 
#' Doesn't yet account for missing data.
#' Unintelligent accommodation of "only k' of k states observed"
#' 
#' @examples
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' tree <- TreeTools::NJTree(dataset)
#' SteelInfo(tree, dataset)
#' @importFrom stats splinefun
#' @export
SteelInfo <- function(tree, dataset) {
  lengths <- CharacterLength(tree, dataset, compress = TRUE)
  allLengths <- 0:max(lengths)
  minLengths <- MinimumLength(dataset, compress = TRUE)
  # nTokens <- lengths(apply(PhyDatToMatrix(dataset), 2, intersect, 0:9))
  
  lengthTab <- table(factor(lengths, levels = allLengths), minLengths)

  pTab <- vapply(colnames(lengthTab), function(steps) {
    nTokens <- as.numeric(steps) + 1
    tokenP <- rep(1 / nTokens, nTokens)
    col <- lengthTab[, steps]
    
    maxSteps <- NTip(tree) - nTokens + 1
    lnP <- active_parsimony_dist(tree, tokenP, maxSteps)
    
    lnPFinite <- lnP[is.finite(lnP)]
    nStep <- as.numeric(names(lnPFinite))
    
    sf <- splinefun(nStep, lnPFinite)
    expect <- parsimony_moments(tree, tokenP)
    lnP0 <- sf(expect$expectation)
    #lnP1 <- lnP[[steps]] - lnP0
    (lnP - lnP0)[rownames(lengthTab)]
  }, double(nrow(lengthTab)))
  
  charScore <- pTab[cbind(as.character(lengths), minLengths)][
    attr(dataset, "index")]
  charMax <- pTab[cbind(as.character(minLengths), minLengths)][
    attr(dataset, "index")]
  normScore <- charScore / charMax
  normScore[charMax == 0] <- 1
  
  
  wSum <- sum(normScore * charMax) / sum(charMax)
  structure(
    wSum,
    byChar = normScore,
    charScore = charScore,
    charMax = charMax
  )
}
