#' Rank taxa by their influence on phylogenetic results
#' 
#' `TaxonInfluence()` ranks taxa according to their influence on the optimal
#' topology.  This follows the approach of 
#' \insertCite{Mariadassou2012sb;textual}{TreeSearch} in repeating tree search
#' whilst leaving each taxon in turn out of the analysis, and measuring
#' the distance of reconstructed trees from the optimal tree obtained when
#' all taxa are included in phylogenetic inference.
#' 
#' As \insertCite{Denton2018ee;textual}{TreeSearch} emphasize, the
#' Robinson&ndash;Foulds distance is unsuitable for this purpose; this function
#' allows the user to specify a preferred tree distance measure, defaulting
#' to the clustering information distance \insertCite{Smith2020}{TreeSearch}.
#' Because optimal parsimony trees are not equiprobable, taxon influence is
#' ranked based on the maximum tree-to-tree distance between optimal trees.
#' 
#' @template datasetParam
#' @param tree Optimal tree or summary tree (of class "phylo") or list of trees
#' (of class "list" or "multiPhylo") against which results should be evaluated.
#' If `NULL`, an optimal tree will be sought using parsimony search with 
#' the parameters provided in \dots.
#' @param Distance Function to calculate tree distance; default:
#' [`ClusteringInfoDistance()`].
#' 
#' @param \dots Parameters for [`MaximizeParsimony()`].
#' Tree search will be conducted using `tree` as a starting tree.
#' 
#' @returns `TaxonInfluence()` returns a named vector, listing the phylogenetic
#' influence of each taxon, measured in the units of the chosen tree distance
#' metric (default = bits).
#' 
#' @references \insertAllCited{}
#' 
#' @template MRS
#' @examples
#' #' # Load data for analysis in R
#' library("TreeTools")
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' 
#' # Small dataset for demonstration purposes
#' dataset <- congreveLamsdellMatrices[[42]][1:8, ]
#' 
#' TaxonInfluence(dataset, ratchIter = 0, startIter = 0, verb = 0)
#'
#' @family tree scoring
#' @importFrom TreeDist ClusteringInfoDistance
#' @importFrom cli cli_h1
#' @export
TaxonInfluence <- function(
    dataset,
    tree = NULL,
    Distance = ClusteringInfoDistance,
    ...
  ) {
  if (is.null(tree)) {
    tree <- MaximizeParsimony(dataset, ...)
  }
  startTree <- MakeTreeBinary(if (inherits(tree, "phylo")) {
    tree
  } else {
    tree[[1]]
  })
  if (!inherits(startTree, "phylo")) {
    stop("`tree` must be an object / list of objects of class \"phylo\"")
  }
  
  # Return:
  vapply(names(dataset), function(leaf) {
    if (is.missing(verbosity) || 
        verbosity > 0) {
      cli_h1(paste("Taxon influence:", leaf))
    }
    result <- MaximizeParsimony(
      dataset = dataset[setdiff(names(dataset), leaf)],
      tree = DropTip(startTree, leaf),
      ...
    )
    max(Distance(tree, result))
  }, double(1))
}
