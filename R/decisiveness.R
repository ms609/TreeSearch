#' Decisivness
#' 
#' Calculate the "Information Decisiveness" statistic of Goloboff (1991)
#' 
#' "Information decisiveness" is the extent to which a dataset supports "a
#' choice or decision between different classifications" (Goloboff 1991).
#' 
#' A dataset with a high decisiveness provives a good basis for preferring
#' certain trees over others.
#' 
#TODO: BUT SEE: De Laet & Smets 1999, Data decisiveness, missing entries, and the DD index
#TODO add to pkgdown index
#' @param 
#' 
#' @return `Decisiveness()` returns the value of Goloboff's (1991) decisiveness
#' statistic, a measure of the extent to which a dataset provides strong support
#' to its associated most parsimonious trees.
#' 
#' @references 
#' \insertRef{Goloboff1991}{TreeSearch}
#' @export
Decisiveness <- function (dataset, bestScore) {
  # sBar = mean number of steps for a character on a tree. See Archie & Felsenstein 1993 for (01) formula.
  
  s <- bestScore
  m <- sum(MinimumSteps(dataset))
  
  # Return:
  (sBar - s) / (sBar - m)
  
}