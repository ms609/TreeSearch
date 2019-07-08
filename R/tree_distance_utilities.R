#' Wrapper for tree distance calculations
#' 
#' Calls tree distance functions from trees or lists of trees
#' 
#' @inheritParams MutualArborealInfo
#' @param Func Tree distance function.
#' @param \dots Additional arguments to `Func`.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
CalculateTreeDistance <- function (Func, tree1, tree2, reportMatching, ...) {
  if (class(tree1) == 'phylo') {
    if (class(tree2) == 'phylo') {
      if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Tree tips must bear identical labels")
      }
      
      Func(Tree2Splits(tree1), Tree2Splits(tree2), reportMatching, ...)
    } else {
      splits1 <- Tree2Splits(tree1)
      vapply(tree2, 
             function (tr2) Func(splits1, Tree2Splits(tr2), ...),
             double(1))
    }
  } else {
    if (class(tree2) == 'phylo') {
      splits1 <- Tree2Splits(tree2)
      vapply(tree1, 
             function (tr2) Func(splits1, Tree2Splits(tr2), ...),
             double(1))
    } else {
      splits1 <- lapply(tree1, Tree2Splits)
      splits2 <- lapply(tree2, Tree2Splits)
      matrix(mapply(Func, rep(splits1, each=length(splits2)), splits2), 
             length(splits2), length(splits1),
             dimnames = list(names(tree2), names(tree1)), ...)
    }
  }
}

#' Normalize information against total present in both starting trees
#' @param unnormalized Numeric value to be normalized.
#' @param tree1,tree2 Trees from which `unnormalized` was calculated
#' @param InfoInTree Function to calculate the information content of each tree
#' @param how Method for normalization
#' @param Func Function that takes as inputs `tree1Info` and `tree2Info`, and
#' returns a normalizing constant against which to divide `unnormalized`.
#' @param \dots Additional parameters to `InfoInTree`` or `how`.
#' @keywords internal
#' @author Martin R. Smith
#' @export
NormalizeInfo <- function (unnormalized, tree1, tree2, InfoInTree, 
                           how = TRUE, Combine = sum, ...) {
  if (mode(how) == 'logical') {
    tree1Info <- InfoInTree(tree1, ...)
    tree2Info <- InfoInTree(tree2, ...)
  } else if (mode(how) == 'function') {
    tree1Info <- how(tree1, ...)
    tree2Info <- how(tree2, ...)
  } else {
    return(unnormalized / how)
  }
  if (length(tree1Info) == 1 || length(tree2Info) == 1) {
    unnormalized / mapply(Combine, tree1Info, tree2Info)
  } else {
    unnormalized / 
      matrix(mapply(Combine, rep(tree2Info, each=length(tree1Info)), tree1Info), 
           length(tree1Info), length(tree2Info))
  }
}