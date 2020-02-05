#' @name AsBinary
#' @aliases AsBinary
#' @title Convert a number to binary
#' @description Provides a (reversed) binary representation of a decimal integer
#' @usage AsBinary(x)
#'
#' @param x Decimal integer to be converted to binary bits
#' 
#' @details 
#' Provides an array corresponding to binary digits 1, 2, 4, 8, 16, ...
#' 
#' Binary number 0100 (= decimal 4) will be represented as 0 0 1.
#' 
#' @return 
#' An array corresponding to binary digits 1, 2, 4, 8, 16, ...
#' 
#' 'Leading zeros' are not included.
#' 
#' @author Martin R. Smith, adapted from code posted to R mailing list by Spencer Graves
#' 
#' @examples
#'   AsBinary(4)  # 0 0 1
#'   AsBinary(10) # 0 1 0 1
#' 
#' @export
AsBinary <- function(x) {
  N <- length(x)
  xMax <- max(x)
  ndigits <- max(1, (floor(logb(xMax, base=2)) + 1L))
  Base.b <- array(NA, dim=c(N, ndigits))
  for (i in 1:ndigits){#i <- 1
    Base.b[, i] <- (x %% 2)
    x <- (x %/% 2)
  }
  if(N == 1) Base.b[1, ] else Base.b
}