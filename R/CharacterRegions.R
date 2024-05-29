#' Character regions
#' 
#' Distinct character regions defined by a tree
#' 
#' @inheritParams PlotCharacter
#' @template MRS
#' @export
CharacterRegions <- function (tree, dataset, method = c("ACCTRAN", "DELTRAN")) {
  
  # Reconcile labels
  datasetTaxa <- names(dataset)
  treeTaxa <- tree[["tip.label"]]
  if(!all(treeTaxa %fin% datasetTaxa)) {
    stop("Taxa in tree missing from dataset:\n  ",
         paste0(setdiff(treeTaxa, datasetTaxa), collapse = ", "))
  }
  dataset <- dataset[treeTaxa]
  
  character_regions(tree, dataset,
                    match.arg(method, c("ACCTRAN", "DELTRAN")) == "ACCTRAN")
    })[attr(dataset, "index")]
}
