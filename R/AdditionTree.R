#' Addition tree
#' 
#' Generates a starting tree by adding each taxon in turn to the most
#' parsimonious location.
#' 
#' @template datasetParam
#' @param sequence Character vector listing sequence in which to add taxa.
#' Randomized if not provided.
#' @template concavityParam
#' @examples 
#' data('Lobo', package = 'TreeTools')
#' AdditionTree(Lobo.phy, concavity = 10)
#' @template MRS
#' @return `AdditionTree()` returns a tree of class `phylo`, rooted on
#' `sequence[1]`.
#' @importFrom TreeTools AddTipEverywhere PectinateTree
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @export
AdditionTree <- function (dataset, sequence, concavity = Inf) {
  taxa <- names(dataset)
  if (missing(sequence)) {
    sequence <- taxa[1]
  }
  nTaxa <- length(taxa)
  if (length(taxa) < 4) {
    return(PectinateTree(taxa))
  }
  unlisted <- setdiff(taxa, sequence)
  if (length(unlisted) > 0) {
    sequence <- c(sequence, sample(unlisted))
  }
  
  tree <- PectinateTree(sequence[1:3])
  cli_progress_bar('Addition tree', total = sum(2 * (4:nTaxa) - 5),
                   auto_terminate = FALSE)
  for (addition in sequence[4:nTaxa]) {
    candidates <- AddTipEverywhere(tree, addition)
    cli_progress_update(length(candidates))
    scores <- TreeLength(candidates, dataset[candidates[[1]]$tip.label], concavity)
    minScore <- which.min(scores)
    nMin <- length(minScore)
    if (nMin > 1) {
      minScore <- minScore[sample.int(nMin, 1)]
    }
    tree <- candidates[[minScore]]
  }
  cli_progress_done()
  
  tree
}