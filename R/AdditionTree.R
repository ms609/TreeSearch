#' Addition tree
#' 
#' Generates a starting tree by adding each taxon in turn to the most
#' parsimonious location.
#' 
#' @inheritParams MaximizeParsimony
#' @param sequence Character or numeric vector listing sequence in which to add
#' taxa. Randomized if not provided.
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' AdditionTree(inapplicable.phyData[["Longrich2010"]], concavity = 10)
#' @template MRS
#' @return `AdditionTree()` returns a tree of class `phylo`, rooted on
#' `sequence[1]`.
#' @importFrom TreeTools AddUnconstrained AddTipEverywhere MatrixToPhyDat
#' PectinateTree
#' @importFrom cli cli_progress_bar cli_progress_update
#' @family tree generation functions
#' @seealso 
#' 
#' Impose a constraint: [`TreeTools::ImposeConstraint()`](
#' https://ms609.github.io/TreeTools/reference/ImposeConstraint)
#' 
#' Neighbour-joining trees: [`TreeTools::NJTree()`](
#' https://ms609.github.io/TreeTools/reference/NJTree.html);
#' [`TreeTools::ConstrainedNJ()`](
#' https://ms609.github.io/TreeTools/reference/ConstrainedNJ)
#' @export
AdditionTree <- function (dataset, concavity = Inf, constraint, sequence) {
  
  # Initialize missing parameters
  taxa <- names(dataset)
  if (missing(sequence)) {
    sequence <- taxa[[1]]
  } else if (is.numeric(sequence)) {
    sequence <- taxa[sequence]
  }
  
  nTaxa <- length(taxa)
  if (length(taxa) < 4) {
    return(PectinateTree(taxa))
  }
  
  unlisted <- setdiff(taxa, sequence)
  if (length(unlisted) > 0) {
    sequence <- c(sequence, sample(unlisted))
  }
  if (!missing(constraint)) {
    constraint <- AddUnconstrained(constraint, taxa)
  }
  
  # PrepareDataXXX attributes only valid for full dataset
  attr(dataset, "info.amounts") <- NULL
  attr(dataset, "min.length") <- NULL
  attr(dataset, "informative") <- NULL
  attr(dataset, "originalIndex") <- NULL
  
  # Starting tree, rooted on first element in sequence
  tree <- PectinateTree(sequence[1:3])
  
  cli_progress_bar("Addition tree", total = sum(2 * (4:nTaxa) - 5))
  for (addition in sequence[4:nTaxa]) {
    candidates <- AddTipEverywhere(tree, addition)
    nCands <- length(candidates)
    
    theseTaxa <- candidates[[1]][["tip.label"]]
    theseData <- .Recompress(dataset[theseTaxa])
    if (is.finite(concavity)) {
      theseData <- PrepareDataIW(theseData)
    } else if (is.character(concavity)) {
      theseData <- suppressMessages(PrepareDataProfile(theseData))
    }
    
    if (!missing(constraint)) {
      if (!inherits(constraint, "phyDat")) {
        if (is.numeric(constraint) && is.null(dim(constraint))) {
          constraint <- t(constraint)
        }
        constraint <- MatrixToPhyDat(t(as.matrix(constraint)))
      }
      thisConstr <- constraint[theseTaxa]
      if (.ConstraintConstrains(thisConstr)) {
        # Constraint constrains theseTaxa
        
        morphyConstr <- PhyDat2Morphy(thisConstr)
        # Calculate constraint minimum score
        constraintLength <- sum(MinimumLength(thisConstr, compress = TRUE) *
                                attr(thisConstr, "weight"))
        
        .Forbidden <- function (edges) {
          preorder_morphy(edges, morphyConstr) != constraintLength
        }
        
      
        candidates <- candidates[!vapply(lapply(candidates, `[[`, "edge"),
                                         .Forbidden, logical(1))]
        UnloadMorphy(morphyConstr)
      }
    }
    
    # Score remaining candidates
    if (length(theseData)) {
      scores <- TreeLength(candidates, theseData, concavity)
      minScore <- which.min(scores)
      nMin <- length(minScore)
      if (nMin > 1) {
        minScore <- minScore[sample.int(nMin, 1)]
      }
      tree <- candidates[[minScore]]
    } else {
      tree <- sample(candidates, 1)[[1]]
    }
    cli_progress_update(nCands)
  }
  tree
}


.ConstraintConstrains <- function(constraint) {
  if (length(constraint[[1]]) < 1) {
    FALSE
  } else {
    contrast <- attr(constraint, "contrast")
    if (dim(contrast)[[2]] < 2) {
      FALSE
    } else {
      cont <- `mode<-`(contrast, "logical")
      nLevel <- dim(contrast)[[1]]
      # Could be > 2× more efficient using lower.tri
      exclude <- vapply(seq_len(nLevel), function(i) {
        colSums(apply(cont, 1, `&`, cont[i, ])) == 0
      }, logical(nLevel))
      
      # TODO Validate; passes existing tests, but these do not include all 
      # edge cases, e.g. 02 03 1 1
      splits <- exclude * tabulate(unlist(constraint), nLevel)
      any(splits[lower.tri(splits)] > 1 & t(splits)[lower.tri(splits)] > 1)
    }
  }
}


.Recompress <- function(dataset) {
  MatrixToPhyDat(PhyDatToMatrix(dataset))
}
