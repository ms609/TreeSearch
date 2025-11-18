#' Add log probabilities
#' 
#' Add log probabilities whilst retaining precision
#' 
#' @param x,y,z Numeric; log-probabilities to add
LogAdd <- function(x, y, z = NULL) {
  if (is.null(z)) {
    ifelse(!is.finite(x), y,
           ifelse(!is.finite(y), x, {
                  maxLog <- pmax(x, y)
                  maxLog + log(exp(x - maxLog) + exp(y - maxLog))
           }))
  } else {
    ifelse(!is.finite(x), LogAdd(y, z),
           ifelse(!is.finite(y), LogAdd(x, z),
                  ifelse(!is.finite(z), LogAdd(x, y), {
                    maxLog <- pmax(x, y, z)
                    maxLog + log(exp(x - maxLog) + exp(y - maxLog) +
                                   exp(z - maxLog))
           })))
  }
}
