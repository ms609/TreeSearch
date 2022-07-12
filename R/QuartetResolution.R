#' Resolution of quartet
#' 
#' @param trees A list of trees of class phylo, or a multiPhylo object
#' @param tips Vector specifying four tips whose relationship should be
#' reported, in a format accepted by \code{\link[TreeTools]{KeepTip}()}.
#' 
#' @return A vector specifying an integer, for each tree, which of `tips` is
#' most closely related to `tips[4]`.
#' 
#' @importFrom TreeTools as.Split KeepTip PolarizeSplits
#' @examples 
#' trees <- inapplicable.trees[["Vinther2008"]]
#' tips <- c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")
#' QuartetResolution(trees, tips)
#' @export
QuartetResolution <- function(trees, tips) {
  fours <- as.integer(vapply(
    lapply(as.Splits(lapply(trees, KeepTip, tips)),
           PolarizeSplits),
    as.raw, raw(1)))
  log2(fours - 1L)
}
