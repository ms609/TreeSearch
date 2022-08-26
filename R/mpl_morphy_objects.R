#' Details the attributes of a morphy object
#'
#' @param object A Morphy object
#' @param \dots any other parameters... 
#'
#' @return A list detailing the number of taxa, internal nodes, and characters and their weights.
#'
#' @template MRS
#' @method summary morphyPtr
#' @family Morphy API functions
#' @importFrom Rcpp compileAttributes
#' @export
summary.morphyPtr <- function (object, ...) {
  ans <- list()
  class(ans) <- "summary.morphyPtr"
  nTax <- mpl_get_numtaxa(object)
  nChar <- mpl_get_num_charac(object)
  charWeights <- MorphyWeights(object)

  ans$nTax <- nTax 
  ans$nChar <- nChar 
  ans$nInternal <- mpl_get_num_internal_nodes(object)
  ans$charWeights <- charWeights
  ans$allStates <- mpl_get_symbols(object)
  # Return:
  ans
}

#' Set and get the character weightings associated with a Morphy object.
#'
#' `MorphyWeights()` details the approximate and exact weights associated with
#' characters in a `Morphy` object; `SetMorphyWeights()` edits them.
#'
#' @template morphyObjParam
#' @return `MorphyWeights()` returns a data frame with two named rows and 
#' one column per character pattern:
#' row 1, `approx`, is a list of integers specifying the approximate (integral)
#' weights used by MorphyLib;
#' row 2, `exact`, is a list of numerics specifying the exact weights specified
#' by the user.
#' 
#' @examples
#' tokens <- matrix(c(
#'   0, 0, 0, 1, 1, 2,
#'   0, 0, 0, 0, 0, 0), byrow = TRUE, nrow = 2L,
#'   dimnames = list(letters[1:2], NULL))
#' pd <- TreeTools::MatrixToPhyDat(tokens)
#' morphyObj <- PhyDat2Morphy(pd)
#' MorphyWeights(morphyObj)
#' if (SetMorphyWeights(c(1, 1.5, 2/3), morphyObj) != 0L) message("Errored")
#' MorphyWeights(morphyObj)
#' morphyObj <- UnloadMorphy(morphyObj)
#' @template MRS
#' @family Morphy API functions
#' @export
MorphyWeights <- function (morphyObj) {
 vapply(seq_len(mpl_get_num_charac(morphyObj)), mpl_get_charac_weight, 
        list('approx' = 0L, 'exact' = 0), morphyobj = morphyObj)
}

#' @rdname MorphyWeights
#' @param weight A vector listing the new weights to be applied to each character
#' @param checkInput Whether to sanity-check input data before applying. 
#'         Defaults to `TRUE` to protect the user from crashes.
#'
#' @return `SetMorphyWeights()` returns the Morphy error code generated when
#' applying `weight`.
#' @export
SetMorphyWeights <- function (weight, morphyObj, checkInput = TRUE) {
  if (checkInput) if (length(weight) != mpl_get_num_charac(morphyObj)) {
    stop("Number of weights not equal to number of character entries.")
  }
  errors <- vapply(seq_along(weight), 
                   function (i) mpl_set_charac_weight(i, weight[i], morphyObj),
                   integer(1))
  if(any(errors != 0)) warning("Morphy Error encountered: ", 
                               mpl_translate_error(errors[errors < 0]))
  mpl_apply_tipdata(morphyObj)
}

#' Read how a Morphy Object handles the inapplicable token
#' 
#' Gaps represented by the inapplicable token can be treated as 'missing data',
#' i.e. as equivalent to the ambiguous token `?`; as an extra state, equivalent
#' to other states such as `0` or `1`; or as 'inapplicable data' using the 
#' algorithm of Brazeau, Guillerme and Smith (2019).
#' 
#' @template morphyObjParam
#' 
#' @return `GapHandler()` returns a character string stating how
#' gaps are handled by `morphyObj`.
#' @examples 
#' morphyObj <- SingleCharMorphy('-0-0', 'Extra')
#' GapHandler(morphyObj)
#' morphyObj <- UnloadMorphy(morphyObj)
#' @family Morphy API functions
#' @template MRS
#' @export
GapHandler <- function (morphyObj) {
  if (!is.morphyPtr(morphyObj)) {
    stop("`morphyObj` is not a valid Morphy object")
  }
    
  ret <- c('Inapplicable', 'Missing data', 'Extra state', 'Unspecified') # 4 = GAP_MAX
  handler <- mpl_get_gaphandl(morphyObj)
  if (handler < 0L) {
    stop("Morphy object error: ", mpl_translate_error(handler))
  }
  # Return:
  ret[1L + handler]
}

#' Initialize a Morphy object from a `phyDat` object
#' 
#' Creates a new Morphy object with the same size and characters as the 
#' `phyDat` object. 
#' Once finished with the object, it should be destroyed using
#' [`UnloadMorphy()`] to free the allocated memory.
#' 
#'
#' @param phy An object of \pkg{phangorn} class \code{phyDat}.
#' @template gapParam
#' @return `PhyDat2Morphy()` returns a pointer to an initialized Morphy object.
#' 
#' @examples
#' data('Lobo', package='TreeTools')
#' morphyObj <- PhyDat2Morphy(Lobo.phy)
#' # Set object to be destroyed at end of session or closure of function
#' # on.exit(morphyObj <- UnloadMorphy(morphyObj), add = TRUE)
#' 
#' # Do something with pointer
#' # ....
#' 
#' # Or, instead of on.exit, manually destroy morphy object and free memory:
#' morphyObj <- UnloadMorphy(morphyObj)
#' @template MRS
#' @family Morphy API functions
#' @importFrom TreeTools PhyToString
#' @export
PhyDat2Morphy <- function (phy, gap = "inapplicable") {
  
  if (!inherits(phy, "phyDat")) {
    stop("Invalid data type ", class(phy), "; should be phyDat.")
  }
  
  morphyObj <- structure(mpl_new_Morphy(), class = "morphyPtr")
  nTax <- length(phy)
  weight <- attr(phy, "weight")
  nChar <- attr(phy, "nr")
  
  if (mpl_init_Morphy(nTax, nChar, morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), " in mpl_init_Morphy")           #nocov
  }
  if (mpl_set_gaphandl(.GapHandler(gap), morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), " in mpl_set_gaphandl")          #nocov
  }
  if (mpl_attach_rawdata(PhyToString(phy, ps=";", useIndex = FALSE,
                                     byTaxon = TRUE, concatenate = TRUE),
                          morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), " in mpl_attach_rawdata")        #nocov
  }
  if (mpl_set_num_internal_nodes(nTax - 1L, morphyObj) -> error) { # One is the "dummy root"
    stop("Error ", mpl_translate_error(error), " in mpl_set_num_internal_nodes")
  }
  if (any(vapply(seq_len(nChar), 
                 function (i) mpl_set_parsim_t(i, "FITCH", morphyObj),
                 NA_integer_) -> error)) {
      stop("Error ", mpl_translate_error(min(error)), "in mpl_set_parsim_t")    #nocov
  }
  if (any(vapply(seq_len(nChar),
                 function (i) mpl_set_charac_weight(i, weight[i], morphyObj),
                 NA_integer_) -> error)) {
    stop("Error ", mpl_translate_error(min(error)), "in mpl_set_charac_weight") #nocov
  }
  if (mpl_apply_tipdata(morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), "in mpl_apply_tipdata")          #nocov
  }
  # Return:
  morphyObj
}

#' Translate a gap treatment into a string in the format expected by Morphy
#' @param gap Character vector: how should gaps be handled?
#' @return Character string that can be translated into a gap handling strategy
#' by Morphy.
#' @keywords internal
.GapHandler <- function (gap) {
  handler <- pmatch(tolower(gap),
                    c("inapplicable", "missing", "ambiguous", "extra state"))
  if (is.na(handler)) {
    stop("`treatment` must be an abbreviation of \"inapplicable\", ",
         "\"missing\" or \"extra state\"")
  }
  
  # Return:
  switch(handler, "inapplicable", "missing", "missing", "newstate")
}

#' Check for error whilst modifying Morphy object
#' @param action action to perform
#'
#' @family Morphy API functions
#' @keywords internal
#' @export
MorphyErrorCheck <- function (action) {
  if (action -> error) {
    stop("Morphy object encountered error ", mpl_translate_error(error), "\n")
  }
}

#' Morphy object from single character
#' 
#' @param char State of each character at each tip in turn, in a format that will be converted
#'             to a character string by \code{\link{paste0}(char, ';', collapse='')}.
#' @template gapParam
#'
#' @return A pointer to an object of class `morphyObj`.
#' Don't forget to unload it when you've finished with it.
#'
#' @examples 
#' morphyObj <- SingleCharMorphy('-0-0', gap = 'Extra')
#' RandomTreeScore(morphyObj)
#' morphyObj <- UnloadMorphy(morphyObj)
#' @template MRS
#' @seealso 
#' Score a tree: [`MorphyTreeLength()`]
#' 
#' @family Morphy API functions
#' @export
SingleCharMorphy <- function (char, gap = "inapp") {
  char <- paste0(c(char, ";"), collapse = "")
  entries <- gregexpr("\\{[^\\{]+\\}|\\([^\\()]+\\)|[^;]", char)
  nTip <- length(entries[[1]])
  morphyObj <- mpl_new_Morphy()
  MorphyErrorCheck(mpl_init_Morphy(nTip, 1, morphyObj)) 
  MorphyErrorCheck(mpl_set_gaphandl(.GapHandler(gap), morphyObj))
  MorphyErrorCheck(mpl_attach_rawdata(char, morphyObj)) 
  MorphyErrorCheck(mpl_set_num_internal_nodes(nTip - 1L, morphyObj)) 
  MorphyErrorCheck(mpl_set_parsim_t(1, "FITCH", morphyObj))
  MorphyErrorCheck(mpl_set_charac_weight(1, 1, morphyObj)) 
  MorphyErrorCheck(mpl_apply_tipdata(morphyObj))
  class(morphyObj) <- "morphyPtr"
  morphyObj
}

#' Is an object a valid Morphy object?
#' @template morphyObjParam
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
#' Destroys a previously-created Morphy object.
#'
#' Best practice is to call `morphyObj <- UnloadMorphy(morphyObj)`
#' Failure to do so will cause a crash if `UnloadMorphy()` is called on an
#' object that  has already been destroyed
#'
#' @template morphyObjParam
#' @return Morphy error code, decipherable using \code{\link{mpl_translate_error}}
#' @author Martin R. Smith
#' @family Morphy API functions
#' @export
UnloadMorphy <- function (morphyObj) {
  if (!is.morphyPtr(morphyObj)) {
    stop ("Object is not a valid pointer; cannot destroy.")
  }
  if (mpl_delete_Morphy(morphyObj) -> error) {
    stop("Error ", mpl_translate_error(error), "in mpl_delete_Morphy")          #nocov
  }
  return (error)
}
