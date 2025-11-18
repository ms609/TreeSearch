#' Add log probabilities
#' 
#' Add log probabilities whilst retaining precision
#' 
#' @param x Numeric vector giving log-probabilities to add
#' @examples
#' exp(LogSumExp(log(1:3)))
#' exp(sum(log(1:3)))
#' 
#' # Intended for large-scale numbers, e.g.
#' choose(2345, 456) # Overflows
#' LogSumExp(lchoose(2345, 456), lchoose(2345, 456))
#' @export
LogSumExp <- function(...) {
  x <- as.numeric(c(...))
  if (all(is.infinite(x) & x == -Inf)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}
