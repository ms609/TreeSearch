#' Character regions
#' 
#' Distinct character regions defined by a tree
#' 
#' @inheritParams PlotCharacter
#' @template MRS
#' @importFrom TreeTools MatrixToPhyDat PhyDatToMatrix
#' @export
CharacterRegions <- function (tree, dataset, method = c("ACCTRAN", "DELTRAN")) {
  
  if (!inherits(dataset, "phyDat")) {
    dataset <- MatrixToPhyDat(dataset)
  }
  
  # Reconcile labels
  datasetTaxa <- names(dataset)
  treeTaxa <- tree[["tip.label"]]
  if(!all(treeTaxa %fin% datasetTaxa)) {
    stop("Taxa in tree missing from dataset:\n  ",
         paste0(setdiff(treeTaxa, datasetTaxa), collapse = ", "))
  }
  dataset <- dataset[treeTaxa]
  
  index <- attr(dataset, "index")
  phyMat <- do.call(rbind, 
                    dataset[, match(seq_len(attr(dataset, "nr")),
                                    attr(dataset, "index"))])
  
  contrast <- attr(dataset, "contrast")
  if (ncol(contrast) > 30) {
    # Won't fit in an int
    stop("Too many columns in contrast matrix")
  }
  binaries <- apply(contrast == 1, 1, function(x) sum(2 ^ (which(x) - 1)))
  inputState <- matrix(binaries[phyMat], nrow(phyMat), ncol(phyMat))
  
  # Return:
  lapply(character_regions(
    tree, inputState, match.arg(method, c("ACCTRAN", "DELTRAN")) == "ACCTRAN"
  ), sort, decreasing = TRUE)[index]
}
