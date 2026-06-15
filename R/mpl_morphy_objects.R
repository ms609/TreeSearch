# Morphy objects in TreeSearch 2.0 are lightweight, MorphyLib-free handles.
# A `morphyPtr` carries a `"native"` attribute (see `.NativeMorphyData()`) that
# holds everything the native Fitch kernel needs to score a tree, including a
# mutable character-weight vector used by the custom-search resampling
# functions (`Jackknife()`, `MorphyBootstrap()`).

#' Initialize a Morphy object from a `phyDat` object
#'
#' Creates a lightweight scoring handle with the same characters as the
#' `phyDat` object, for use with the custom-search functions ([`TreeSearch()`],
#' [`Ratchet()`], [`Jackknife()`]).
#' Once finished with the object, it should be released using [`UnloadMorphy()`].
#'
#' @param phy An object of \pkg{phangorn} class \code{phyDat}.
#' @param gap Character specifying how gaps (the inapplicable token) will be
#' handled.  Only the default `"inapplicable"` treatment (the algorithm of
#' Brazeau, Guillerme and Smith 2019) is supported; recode your data if you
#' require gaps to be treated as missing or as an extra state.
#' @param weight Numeric giving weight to apply to each character.
#' Will be recycled.
#' @return `PhyDat2Morphy()` returns an object of class `morphyPtr`.
#'
#' @examples
#' data("Lobo", package="TreeTools")
#' morphyObj <- PhyDat2Morphy(Lobo.phy)
#'
#' # Do something with the object
#' # ....
#'
#' # Release the object when finished
#' morphyObj <- UnloadMorphy(morphyObj)
#' @template MRS
#' @family Morphy API functions
#' @useDynLib TreeSearch, .registration = TRUE
#' @importFrom Rcpp compileAttributes
#' @export
PhyDat2Morphy <- function(phy, gap = "inapplicable",
                          weight = attr(phy, "weight")) {

  if (!inherits(phy, "phyDat")) {
    stop("Invalid data type ", class(phy), "; should be phyDat.")
  }
  native <- .NativeMorphyData(phy, gap, weight)
  if (is.null(native)) {
    if (.GapHandler(gap) != "inapplicable") {
      stop("Only the `inapplicable` gap treatment is supported; ",
           "recode your data to treat gaps as missing or as an extra state.")
    }
    stop("Could not extract scoring data from `phy`.")
  }
  morphyObj <- structure(
    list(nTip = native[["nTip"]], nChar = length(native[["weight"]])),
    class = "morphyPtr")
  attr(morphyObj, "native") <- native
  # Return:
  morphyObj
}

#' Translate a gap treatment into a canonical string
#' @param gap Character vector: how should gaps be handled?
#' @return Character string identifying the gap-handling strategy.
#' @keywords internal
.GapHandler <- function(gap) {
  handler <- pmatch(tolower(gap),
                    c("inapplicable", "missing", "ambiguous", "extra state"))
  if (is.na(handler)) {
    stop("`treatment` must be an abbreviation of \"inapplicable\", ",
         "\"missing\" or \"extra state\"")
  }

  # Return:
  switch(handler, "inapplicable", "missing", "missing", "newstate")
}

# Precompute MorphyLib-free Fitch scoring data for a Morphy object.
# Returns NULL (no native support) for the rarely used `ambiguous`/`extra` gap
# modes; the default `inapplicable` mode is scored by the native kernel, which
# handles ambiguous-with-inapplicable tokens correctly.
.NativeMorphyData <- function(phy, gap = "inapplicable",
                              weight = attr(phy, "weight")) {
  if (.GapHandler(gap) != "inapplicable") {
    return(NULL)
  }
  at <- attributes(phy)
  if (is.null(at[["levels"]]) || ncol(at[["contrast"]]) == 0L) {
    return(NULL)
  }
  list(
    phy = phy,
    weight = rep_len(weight, at[["nr"]]),
    tip.label = names(phy),
    nTip = length(phy)
  )
}

# As `.NativeMorphyData()`, but reconstructing a single-character phyDat from a
# Morphy token string (e.g. "-0-{12}0") for `SingleCharMorphy()`.
#' @importFrom TreeTools MatrixToPhyDat
.SingleCharNative <- function(char, gap = "inapp") {
  if (.GapHandler(gap) != "inapplicable") {
    return(NULL)
  }
  charStr <- paste0(paste0(char, collapse = ""), ";")
  entries <- regmatches(
    charStr, gregexpr("\\{[^{}]+\\}|\\([^()]+\\)|[^;]", charStr))[[1]]
  nTip <- length(entries)
  if (nTip < 1L) {
    return(NULL)
  }
  m <- matrix(entries, ncol = 1L,
              dimnames = list(as.character(seq_len(nTip)), NULL))
  phy <- tryCatch(MatrixToPhyDat(m), error = function(e) NULL)
  if (is.null(phy)) {
    return(NULL)
  }
  .NativeMorphyData(phy, "inapplicable", attr(phy, "weight"))
}

# Read the character-weight vector from a Morphy object's native data.
.MorphyWeight <- function(morphyObj) {
  native <- attr(morphyObj, "native", exact = TRUE)
  if (is.null(native)) {
    stop("`morphyObj` lacks native scoring data")
  }
  native[["weight"]]
}

# Return a copy of `morphyObj` whose characters carry the supplied weights.
# Used by the custom-search resampling functions to score a resampled dataset
# without mutating the original object.
.SetMorphyWeight <- function(morphyObj, weight) {
  native <- attr(morphyObj, "native", exact = TRUE)
  if (is.null(native)) {
    stop("`morphyObj` lacks native scoring data")
  }
  native[["weight"]] <- rep_len(weight, length(native[["weight"]]))
  attr(morphyObj, "native") <- native
  morphyObj
}

#' Morphy object from single character
#'
#' @param char State of each character at each tip in turn, in a format that
#' will be converted to a character string by
#' \code{\link{paste0}(char, ";", collapse="")}.
#' @inheritParams PhyDat2Morphy
#'
#' @return A `morphyPtr` object.
#' Release it with [`UnloadMorphy()`] when finished.
#'
#' @examples
#' morphyObj <- SingleCharMorphy("-0-0")
#' RandomTreeScore(morphyObj)
#' morphyObj <- UnloadMorphy(morphyObj)
#' @template MRS
#' @seealso
#' Score a tree: [`MorphyTreeLength()`]
#'
#' @family Morphy API functions
#' @export
SingleCharMorphy <- function (char, gap = "inapp") {
  native <- .SingleCharNative(char, gap)
  if (is.null(native)) {
    if (.GapHandler(gap) != "inapplicable") {
      stop("Only the `inapplicable` gap treatment is supported; ",
           "recode your data to treat gaps as missing or as an extra state.")
    }
    stop("Could not parse character `", paste0(char, collapse = ""), "`.")
  }
  morphyObj <- structure(
    list(nTip = native[["nTip"]], nChar = length(native[["weight"]])),
    class = "morphyPtr")
  attr(morphyObj, "native") <- native
  # Return:
  morphyObj
}

#' Is an object a valid Morphy object?
#' @inheritParams MorphyTreeLength
#' @return `is.morphyPtr()` returns `TRUE` if `morphyObj` is a valid morphy
#' pointer, `FALSE` otherwise.
#' @template MRS
#' @family Morphy API functions
#' @export
is.morphyPtr <- function (morphyObj) {
  inherits(morphyObj, "morphyPtr")
}

#' Destroy a Morphy object
#'
#' Releases a previously-created Morphy object.
#' Retained for backwards compatibility: Morphy objects in TreeSearch 2.0 are
#' managed by R's garbage collector, so `UnloadMorphy()` simply invalidates the
#' handle and returns `0`.
#'
#' @inheritParams MorphyTreeLength
#' @return `UnloadMorphy()` returns the integer `0`, invisibly.
#' @author Martin R. Smith
#' @family Morphy API functions
#' @export
UnloadMorphy <- function (morphyObj) {
  if (!is.morphyPtr(morphyObj)) {
    stop ("Object is not a valid pointer; cannot destroy.")
  }
  # Return:
  invisible(0L)
}
