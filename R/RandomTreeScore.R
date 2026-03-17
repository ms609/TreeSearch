#' Parsimony score of random tree
#' 
#' Generates a random tree topology and returns its parsimony score under
#' equal weights.
#' 
#' @param dataset A phyDat object (recommended) or a Morphy object created
#'   with [`PhyDat2Morphy()`] (legacy; deprecated).
#'
#' @return `RandomTreeScore()` returns a numeric parsimony score.
#' @examples 
#' tokens <- matrix(c(
#'   0, "-", "-", 1, 1, 2,
#'   0, 1, 0, 1, 2, 2,
#'   0, "-", "-", 0, 0, 0), byrow = TRUE, nrow = 3L,
#'   dimnames = list(letters[1:3], NULL))
#' pd <- TreeTools::MatrixToPhyDat(tokens)
#' RandomTreeScore(pd)
#' @importFrom TreeTools RandomTree
#' @export
RandomTreeScore <- function(dataset) {
  if (inherits(dataset, "morphyPtr")) {
    nTip <- mpl_get_numtaxa(dataset)
    if (nTip < 2) {
      return(0L)
    }
    return(.Call(`RANDOM_TREE_SCORE`, as.integer(nTip), dataset))
  }
  
  nTip <- length(dataset)
  if (nTip < 2) {
    return(0)
  }
  tree <- RandomTree(dataset, root = TRUE)
  TreeLength(tree, dataset)
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
