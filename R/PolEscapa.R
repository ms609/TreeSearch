#' Contribution of character to leaf instability
#' 
#' Would tree length change if a character was coded as ambiguous for each 
#' leaf \insertRef{Pol2009}{TreeSearch}?
#' 
#' @param trees List of trees of class `phylo`, or `multiPhylo` object.
#' @param char `phyDat` object containing a single character.
#' @template concavityParam
#' 
#' @return `PolEscapa()` returns a named numeric vector listing the mean change
#' to tree length resulting if the character were coded ambiguous for each
#' leaf in turn, under the specified concavity constant.
#' 
#' @references \insertAllCited{}
#' @template MRS
#' @examples
#' trees <- inapplicable.trees[["Vinther2008"]]
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' char <- dataset[, 10]
#' PolEscapa(trees, dataset[, 10])
#' @export
PolEscapa <- function(trees, char, concavity = Inf) {
  if(!inherits(char, "phyDat")) {
    stop("`char` must be a character of class `phyDat`.")
  }
  start <- TreeLength(trees, char, concavity)
  cont <- attr(char, "contrast")
  qm <- which(rowSums(cont) == dim(cont)[2])
  loss <- colSums(
    abs(start - vapply(seq_along(char), function(leaf) {
      charQm <- char
      charQm[[leaf]] <- qm
      TreeLength(trees, charQm, concavity)
    }, start))
  )
  
  # Return:
  loss / length(trees)
}