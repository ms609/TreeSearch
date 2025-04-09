#' Consistency / retention "indices"
#' 
#' `Consistency()` calculates the consistency "index" and retention index
#' \insertCite{Farris1989}{TreeSearch}
#' for each character in a dataset, given a bifurcating tree.
#' Although there is not a straightforward interpretation of these indices,
#' they are sometimes taken as an indicator of the fit of a character to a 
#' tree.
#' Values correlate with the number of species sampled and the
#' distribution of taxa between character states, so are not strictly comparable
#' between characters in which these factors differ.
#' 
#' The **consistency "index"** \insertCite{Kluge1969}{TreeSearch} is defined as the
#' number of steps observed in the most parsimonious mapping of a character
#' to a tree, divided by the number of steps observed on the shortest possible
#' tree for that character. A value of one indicates that a character's fit to
#' the tree is optimal.
#' Note that as the possible values of the consistency index do not range from
#' zero to one, it is not an index in the mathematical sense of the term.
#' 
#' The maximum length of a character (see [`MaximumLength()`]) is the
#' number of steps in a parsimonious reconstruction on the longest possible tree
#' for a character. 
#' The **retention index** is the maximum length of a character minus the number
#' of steps observed on a given tree; divided by the maximum length minus the
#' minimum length.  It is interpreted as the ratio between the observed 
#' homoplasy, and the maximum observed homoplasy, and scales from zero
#' (worst fit that can be reconstructed under parsimony) to one (perfect fit).
#' 
#' The **rescaled consistency index** is the product of the consistency and
#' retention indices; it rescales the consistency index such that its range of
#' possible values runs from zero (least consistent) to one
#' (perfectly consistent).
#' 
#' The lengths of characters including inapplicable tokens are calculated
#' following \insertCite{Brazeau2019;textual}{TreeSearch}, matching their
#' default treatment in [`TreeLength()`].
#' 
#' @inheritParams CharacterLength
#' 
#' @return `Consistency()` returns a matrix with named columns specifying the 
#' consistency index (`ci`),
#' retention index (`ri`), and
#' rescaled consistency index (`rc`).
#' 
#' @examples 
#' data(inapplicable.datasets)
#' dataset <- inapplicable.phyData[[4]]
#' head(Consistency(dataset, TreeTools::NJTree(dataset)))
#' @references \insertAllCited{}
#' @template MRS
#' @export
Consistency <- function (dataset, tree, compress = FALSE) {
  dsTips <- TipLabels(dataset)
  trTips <- TipLabels(tree)
  if (!setequal(dsTips, trTips)) {
    dsHas <- setdiff(dsTips, trTips)
    trHas <- setdiff(trTips, dsTips)
    stop("Tip label mismatch: ",
         if (length(dsHas)) "\n   `dataset` has ", paste(dsHas, collapse = ", "),
         if (length(trHas)) "\n   `tree` has ", paste(trHas, collapse = ", ")
    )
  }
  minLength <- MinimumLength(dataset, compress = TRUE) # Farris's m
  maxLength <- MaximumLength(dataset, compress = TRUE) # Farris's g
  obsLength <- CharacterLength(tree, dataset, compress = TRUE) # farris's s
  extra <- obsLength - minLength # Farris's h
  maxHomoplasy <- (maxLength - minLength) # g - m
  
  ci <- minLength / obsLength # Farris's c = m / s
  distortion <- extra / maxHomoplasy # Farris's d = h / (g - m)
  
  ri <- (maxLength - obsLength) / maxHomoplasy
  
  rc <- ri * minLength / obsLength

  ret <- cbind(ci = ci, ri = ri, rc = rc)
  
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(dataset, "index"), ]
  }
}
