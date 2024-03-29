#' Relationship between four taxa
#' 
#' @param trees A list of trees of class `phylo`, or a `multiPhylo` object.
#' @param tips Vector specifying four tips whose relationship should be
#' reported, in a format accepted by \code{\link[TreeTools]{KeepTip}()}.
#' 
#' @return A vector specifying an integer, for each tree, which of `tips[-1]`
#' is most closely related to `tips[1]`.
#' 
#' @examples 
#' trees <- inapplicable.trees[["Vinther2008"]]
#' tips <- c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")
#' QuartetResolution(trees, tips)
#' @importFrom TreeTools as.Splits KeepTip PolarizeSplits
#' @family utility functions
#' @export
QuartetResolution <- function(trees, tips) {
  fours <- as.integer(vapply(
    lapply(as.Splits(KeepTip(trees, tips), tips), PolarizeSplits),
    as.raw,
    raw(1)
  ))
  log2(fours - 1L)
}
