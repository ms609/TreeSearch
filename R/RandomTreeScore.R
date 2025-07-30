#' Parsimony score of random postorder tree
#' 
#' @inheritParams MorphyTreeLength
#'
#' @return `RandomTreeScore()` returns the parsimony score of a random tree
#'  for the given Morphy object.
#' @examples 
#' tokens <- matrix(c(
#'   0, "-", "-", 1, 1, 2,
#'   0, 1, 0, 1, 2, 2,
#'   0, "-", "-", 0, 0, 0), byrow = TRUE, nrow = 3L,
#'   dimnames = list(letters[1:3], NULL))
#' pd <- TreeTools::MatrixToPhyDat(tokens)
#' morphyObj <- PhyDat2Morphy(pd)
#'
#' RandomTreeScore(morphyObj)
#' 
#' morphyObj <- UnloadMorphy(morphyObj)
#' @export
RandomTreeScore <- function (morphyObj) {
  nTip <- mpl_get_numtaxa(morphyObj)
  if (nTip < 2) {
    # Return:
    0L
  } else {
    # Return:
    .Call(`RANDOM_TREE_SCORE`, as.integer(nTip), morphyObj)
  }
}

#' Random postorder tree
#' 
#' @param nTip Integer specifying the number of tips to include in the tree
#' (minimum 2).
#'
#' @return A list with three elements, each a vector of integers, respectively 
#' containing:
#' 
#'  - The parent of each tip and node, in order
#'         
#'  - The left child of each node
#'         
#'  - The right child of each node.
#'
#' @family tree generation functions
#' @export
RandomMorphyTree <- function (nTip) {  
  if (nTip < 2) {
    stop("nTip < 2 not implemented: a tip is not a tree.")
  }
  # Return:
  .Call(`RANDOM_TREE`, as.integer(nTip))
}
