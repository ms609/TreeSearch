#' When was a tree topology first hit?
#' 
#' Reports when each tree in a list was first found by tree search.
#' This information is read from the `firstHit` attribute if present.
#' If not, the names of the trees are assumed to correspond to the iteration
#' in which they were first hit - the situation when trees found by TreeSearch
#' are saved to file.
#' 
#' @param trees A list of trees, or a `multiPhylo` object.
#' @return `trees`, with a `firstHit` attribute listing the search iteration in
#' which each topology was first found by tree search.
#' @template MRS
#' 
#' @examples
#' treeFile <- system.file("datasets/Wills2012.nex", package = "TreeSearch")
#' trees <- read.nexus(treeFile)
#' names(trees)
#' attr(WhenFirstHit(trees), "firstHit")
#' @export
WhenFirstHit <- function(trees) {
  if (is.null(attr(trees, "firstHit"))) {
    treeNames <- names(trees)
    pattern <- "(seed|start|ratch\\d+|final)_\\d+"
    if (length(grep(pattern, treeNames, perl = TRUE)) == length(trees)) {
      whenHit <- gsub(pattern, "\\1", treeNames, perl = TRUE)
      
      attr(trees, "firstHit") <- table(whenHit)[unique(whenHit)]
    }
  }
  
  # Return:
  trees
}
