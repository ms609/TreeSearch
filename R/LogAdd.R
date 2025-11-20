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
  ret <- m + log(sum(exp(x - m)))
  if (abs(ret) < sqrt(.Machine$double.eps)) 0 else ret
}

#' Multiply log probabilities
#' @param x List of expressions to log-multiply. If any expressions evaluate to
#' `-Inf`, `log(0)` will be returned, without evaluating further expressions.
#' @export
LogProdExp <- function(x) {
  result <- 0
  for (i in seq_along(x)) {
    step <- x[[i]]
    if (step == -Inf) {
      return(-Inf)
    }
    result <- result + step
  }
  result
}
