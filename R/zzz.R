# Suppress "NOTE: Nothing imported from Rdpack":
#' @useDynLib TreeSearch, .registration = TRUE
.onUnload <- function (libpath) {
  library.dynam.unload("TreeSearch", libpath)
}
