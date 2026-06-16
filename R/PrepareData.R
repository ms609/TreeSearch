# TreeSearch 2.0 prepares a dataset for fast, repeated native-C++ parsimony
# scoring as a lightweight `ParsimonyData` object: pre-extracted data matrices
# plus a mutable integer weight vector (used by the custom-search resampling
# functions `Jackknife()` and `BootstrapTree()`).  Equal weights, implied
# weights and profile parsimony are selected through the `concavity` argument.

#' Prepare a dataset for tree scoring
#'
#' `PrepareData()` extracts the matrices that the native C++ Fitch engine needs
#' in order to score a tree repeatedly, for use with the custom-search
#' functions [`TreeSearch()`], [`Ratchet()`] and [`Jackknife()`] (as their
#' `InitializeData` argument).  Once finished, release it with [`ReleaseData()`]
#' (a formality: the object holds no external resources).
#'
#' @param dataset An object of \pkg{phangorn} class `phyDat`.
#' @param concavity Determines the optimality criterion under which trees are
#' scored:
#'  `Inf` (the default) for equal-weights parsimony;
#'  a finite, positive number for implied weighting with that concavity constant
#'  \insertCite{Goloboff1993}{TreeSearch};
#'  or the string `"profile"` for profile parsimony
#'  \insertCite{Faith2001}{TreeSearch}.
#' @return `PrepareData()` returns an object of class `ParsimonyData`: a list of
#'  pre-extracted scoring matrices suitable for [`EdgeListScore()`].
#' @template MRS
#' @seealso
#' Custom searches that consume the object: [`TreeSearch()`], [`Ratchet()`],
#' [`Jackknife()`].
#' @references \insertAllCited{}
#' @family custom search functions
#' @useDynLib TreeSearch, .registration = TRUE
#' @importFrom Rcpp compileAttributes
#' @export
PrepareData <- function(dataset, concavity = Inf) {
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }

  useProfile <- identical(concavity, "profile")
  if (useProfile) {
    dataset <- PrepareDataProfile(dataset)
    infoAmounts <- attr(dataset, "info.amounts")
  } else {
    infoAmounts <- NULL
  }

  iw <- !useProfile && is.finite(concavity)
  if (iw) {
    if (concavity <= 0) {
      stop("`concavity` must be positive, `Inf` for equal weights, ",
           "or \"profile\" for profile parsimony.")
    }
    if (!("min.length" %in% names(attributes(dataset)))) {
      dataset <- PrepareDataIW(dataset)
    }
    minSteps <- as.integer(attr(dataset, "min.length"))
  } else {
    minSteps <- integer(0)
  }

  at <- attributes(dataset)
  # `original_weight` (integer pattern multiplicities) is the base for
  # character resampling; `weight` is the rescaled score weight passed to the
  # C++ engine (identical for the usual integer weights).
  structure(
    list(
      contrast = at[["contrast"]],
      tip_data = matrix(unlist(dataset, use.names = FALSE),
                        nrow = length(dataset), byrow = TRUE),
      weight = .ScaleWeight(at[["weight"]]),
      levels = at[["levels"]],
      min_steps = minSteps,
      concavity = if (iw) as.double(concavity) else Inf,
      info_amounts = infoAmounts,
      original_weight = as.integer(at[["weight"]]),
      index = at[["index"]],
      tip.label = names(dataset),
      nTip = length(dataset)
    ),
    class = "ParsimonyData")
}

#' @describeIn PrepareData Release a prepared dataset.  A no-op retained for the
#'  custom-search `CleanUpData` hook (the object is managed by R's garbage
#'  collector).
#' @param x A `ParsimonyData` object.
#' @return `ReleaseData()` returns its argument invisibly.
#' @export
ReleaseData <- function(x) {
  invisible(x)
}

#' @describeIn PrepareData Test whether an object is a `ParsimonyData` handle.
#' @return `is.ParsimonyData()` returns `TRUE` or `FALSE`.
#' @export
is.ParsimonyData <- function(x) {
  inherits(x, "ParsimonyData")
}

#' Prepare a single character for scoring
#'
#' Builds a one-character [`ParsimonyData`][PrepareData] object, chiefly for
#' [`RandomTreeScore()`] and testing.
#'
#' @param char State of each tip in turn, in a format understood by
#'  [`TreeTools::MatrixToPhyDat()`] (e.g. the Morphy token string `"-0-0"`).
#' @return `SingleCharData()` returns a `ParsimonyData` object.
#' @examples
#' obj <- SingleCharData("-0-0")
#' RandomTreeScore(obj)
#' @template MRS
#' @family custom search functions
#' @importFrom TreeTools MatrixToPhyDat
#' @export
SingleCharData <- function(char) {
  charStr <- paste0(paste0(char, collapse = ""), ";")
  entries <- regmatches(
    charStr, gregexpr("\\{[^{}]+\\}|\\([^()]+\\)|[^;]", charStr))[[1]]
  nTip <- length(entries)
  if (nTip < 1L) {
    stop("Could not parse character `", paste0(char, collapse = ""), "`.")
  }
  m <- matrix(entries, ncol = 1L,
              dimnames = list(as.character(seq_len(nTip)), NULL))
  PrepareData(MatrixToPhyDat(m))
}

# Canonicalize a gap treatment.  Only "inapplicable" is supported natively;
# used by the deprecated `PhyDat2Morphy()`/`SingleCharMorphy()` shims to reject
# the MorphyLib-era missing/extra-state treatments.
.GapHandler <- function(gap) {
  handler <- pmatch(tolower(gap),
                    c("inapplicable", "missing", "ambiguous", "extra state"))
  if (is.na(handler)) {
    stop("`gap` must be an abbreviation of \"inapplicable\", \"missing\" ",
         "or \"extra state\"")
  }
  switch(handler, "inapplicable", "missing", "missing", "newstate")
}
