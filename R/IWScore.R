#' @template pointlessDots
#' @rdname TreeLength
#' @export
IWScore <- function (tree, dataset, concavity = 10L, ...) {
  .Deprecated("TreeLength")
  TreeLength(tree, dataset, concavity)
}

#' @rdname TreeSearch
#' @export
IWTreeSearch <- function (...) {
  .Deprecated("MaximizeParsimony") # Retained as template, for now.
  message("See also the vignette \"custom optimality criteria\"")
}
