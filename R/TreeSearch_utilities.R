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

# Helper function to install suggested packages with user prompt
.InstallSuggestedPackage <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (interactive()) {
      message("Package '", package, "' is required but not installed.")
      response <- readline(prompt = paste0("Install ", package, "? (y/n): "))
      if (tolower(trimws(response)) == "y") {
        install.packages(package)
        if (!requireNamespace(package, quietly = TRUE)) {
          stop("Failed to install package '", package, "'", call. = FALSE)
        }
      } else {
        stop("Package '", package, "' is required. ",
             "Install it with: install.packages('", package, "')",
             call. = FALSE)
      }
    } else {
      stop("Package '", package, "' is required. ",
           "Install it with: install.packages('", package, "')",
           call. = FALSE)
    }
  }
}

