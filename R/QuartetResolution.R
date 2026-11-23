#' Relationship between four taxa
#' 
#' @param trees A list of trees of class `phylo`, or a `multiPhylo` object.
#' @param tips Vector specifying four tips whose relationship should be
#' reported, in a format accepted by \code{\link[TreeTools]{KeepTip}()}.
#' 
#' @return A vector specifying an integer, for each tree, which of `tips[-1]`
#' is most closely related to `tips[1]`.  A tree in which the four tips form
#' an unresolved (star) quartet contributes `NA` to this vector.
#'
#' @examples
#' trees <- inapplicable.trees[["Vinther2008"]]
#' tips <- c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")
#' QuartetResolution(trees, tips)
#' @importFrom TreeTools as.Splits KeepTip PolarizeSplits
#' @family utility functions
#' @export
QuartetResolution <- function(trees, tips) {
  splits <- lapply(as.Splits(KeepTip(trees, tips), tips), PolarizeSplits)
  fours <- unname(vapply(splits, function(x) {
    if (length(x) == 0) {
      NA_integer_ # Unresolved (star) quartet: no split to report
    } else {
      as.integer(as.raw(x))
    }
  }, integer(1)))
  log2(fours - 1L)
}
