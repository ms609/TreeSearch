#' @rdname TreeSearch
#' @return `EmptyPhyDat()` returns a `phyDat` object comprising a single
#' null character, coded with state zero for every leaf in `tree`.
#' @importFrom TreeTools MatrixToPhyDat NTip TipLabels
#' @export
EmptyPhyDat <- function(tree) {
  mat <- matrix(0, NTip(tree), 1, dimnames = list(TipLabels(tree), NULL))
  MatrixToPhyDat(mat)
}

#' @rdname TreeSearch
#' @export
DoNothing <- function(x = NULL) {x}
