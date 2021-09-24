# Suppress "NOTE: Nothing imported from Rdpack":
#' @useDynLib TreeSearch, .registration = TRUE
#' @importFrom TreeTools Preorder
.onUnload <- function (libpath) {
  library.dynam.unload("TreeSearch", libpath)
}
