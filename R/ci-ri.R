#' Consistency / retention 'indices'
#' 
#' `Consistency()` calculates the so-called consistency and retention 'indices'
#' for each character in a dataset, given a bifurcating tree.
#' Although there is not a straightforward interpretation of these indices,
#' they are sometimes taken as an indicator of the fit of a character to a 
#' tree.  Values correlate with the number of species sampled and the 
#' distribution of taxa between character states, so are not strictly comparable
#' between characters in which these factors differ.
#' 
#' #TODO: Retention index not yet implemented.
#' 
#' @template datasetParam
#' @template treeParam
#' 
#' @return `Consistency()` returns a named vector specifying the consistency index (ci),
#' retention index (ri), and rescalced consistency index (rc).
#' 
#' @examples 
#' data(inapplicable.datasets)
#' dataset <- inapplicable.phyData[[4]]
#' Consistency(dataset, TreeTools::NJTree(dataset))
#' @template MRS
#' @export
Consistency <- function (dataset, tree) {
  ci <- MinimumLength(dataset) / CharacterLength(tree, dataset)
  # Return:
  ci
}