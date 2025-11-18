#' Add log probabilities
#' 
#' Add log probabilities whilst retaining precision
#' 
#' @param x,y,z Numeric; log-probabilities to add
#' @export
LogSumExp <- function(x, y, z = NULL) {
  if (is.null(z)) {
    ifelse(!is.finite(x), y,
           ifelse(!is.finite(y), x, {
                  maxLog <- pmax(x, y)
                  maxLog + log(exp(x - maxLog) + exp(y - maxLog))
           }))
  } else if (!is.null(y)) {
    ifelse(!is.finite(x), LogSumExp(y, z),
           ifelse(!is.finite(y), LogSumExp(x, z),
                  ifelse(!is.finite(z), LogSumExp(x, y), {
                    maxLog <- pmax(x, y, z)
                    maxLog + log(exp(x - maxLog) + exp(y - maxLog) +
                                   exp(z - maxLog))
           })))
  } else if (length(x) > 1) {
    LogSumExp(x, LogSumExp(x[-1]))
  } else {
    x
  }
}
