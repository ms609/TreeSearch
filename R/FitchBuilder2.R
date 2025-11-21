#' Count the number of ways to populate skeleton regions with leaves
#'
#' @param regionCounts integer matrix: rows = labels, cols = skeletons.
#'   Each entry = number of regions of that label in the skeleton.
#' @param tokenCount integer vector: number of leaves per label
#' @return integer vector: number of populated trees per skeleton
FitchPopulateCount <- function(regionCounts, tokenCount) {
  
  # Number of rooted trees with n labelled tips
  NRooted <- function(tips) {
    if (tips <= 1) return(1L)
    prod(seq(1, 2 * tips - 3, by = 2))
  }
  
  # Compositions of n into e positive integers
  CompositionsPos <- function(n, e) {
    if (e == 1) return(list(c(n)))
    res <- list()
    k <- 1L
    for (first in 1:(n - e + 1)) {
      tails <- CompositionsPos(n - first, e - 1)
      for (t in tails) {
        res[[k]] <- c(first, t)
        k <- k + 1L
      }
    }
    res
  }
  
  # Count ways to assign leaves to a label with r regions
  LeafAssignment <- function(r, n) {
    if (n < r) return(0L)
    comps <- CompositionsPos(n, r)
    total <- 0
    for (kvec in comps) {
      multinom <- factorial(n) / prod(factorial(kvec))
      total <- total + multinom * prod(vapply(kvec, NRooted, double(1)))
    }
    total
  }
  
  # Apply across skeletons
  apply(regionCounts, 2, function(col) {
    prod(mapply(LeafAssignment, r = col, n = tokenCount))
  })
}

BuildCounter <- function(x, y) {
  FitchBuilder(x, y) |> FitchPopulateCount(y)
}