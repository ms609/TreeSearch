#' Consistency / retention "indices"
#' 
#' `Consistency()` calculates the so-called consistency and retention "indices"
#' \insertCite{Farris1989}{TreeSearch}
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
#' @template compressParam
#' 
#' @return `Consistency()` returns a named matrix with rows specifying the 
#' consistency index (`ci`),
#' retention index (`ri`), and
#' rescaled consistency index (`rc`).
#' 
#' @examples 
#' data(inapplicable.datasets)
#' dataset <- inapplicable.phyData[[4]]
#' Consistency(dataset, TreeTools::NJTree(dataset))
#' @references \insertAllCited{}
#' @template MRS
#' @export
Consistency <- function (dataset, tree, compress = FALSE) {
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
