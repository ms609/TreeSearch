#' Consistency / retention "indices"
#' 
#' `Consistency()` calculates the consistency "index" and retention index
#' \insertCite{Farris1989}{TreeSearch}
#' for each character in a dataset, given a bifurcating tree.
#' Although there is not a straightforward interpretation of these indices,
#' they are sometimes taken as an indicator of the fit of a character to a 
#' tree.
#' Values correlate with the number of species sampled and the
#' distribution of taxa between character states, so are not strictly comparable
#' between characters in which these factors differ; and values cannot be
#' compared between datasets \insertCite{Speed2017}{TreeSearch}.
#' 
#' The **consistency "index"** \insertCite{Kluge1969}{TreeSearch} is defined as the
#' number of steps observed in the most parsimonious mapping of a character
#' to a tree, divided by the number of steps observed on the shortest possible
#' tree for that character. A value of one indicates that a character's fit to
#' the tree is optimal.
#' Note that as the possible values of the consistency index do not range from
#' zero to one, it is not an index in the mathematical sense of the term.
#' Shortcomings of this measure are widely documented 
#' \insertCite{Archie1989,Brooks1986,Steell2023}{TreeSearch}.
#' 
#' The maximum length of a character (see [`MaximumLength()`]) is the
#' number of steps in a parsimonious reconstruction on the longest possible tree
#' for a character. 
#' The **retention index** is the maximum length of a character minus the number
#' of steps observed on a given tree; divided by the maximum length minus the
#' minimum length.  It is interpreted as the ratio between the observed 
#' homoplasy, and the maximum observed homoplasy, and scales from zero
#' (worst fit that can be reconstructed under parsimony) to one (perfect fit).
#' 
#' The **rescaled consistency index** is the product of the consistency and
#' retention indices; it rescales the consistency index such that its range of
#' possible values runs from zero (least consistent) to one
#' (perfectly consistent).
#' 
#' The **relative homoplasy index** \insertCite{Steell2023}{TreeSearch} is
#' the ratio of the observed excess tree length to the excess tree length
#' due to chance, taken as the median score of a character when the leaves
#' of the given tree are randomly shuffled.
#' 
#' The lengths of characters including inapplicable tokens are calculated
#' following \insertCite{Brazeau2019;textual}{TreeSearch}, matching their
#' default treatment in [`TreeLength()`].
#' 
#' @param nRelabel Integer specifying how many times to relabel leaves when
#' estimating null tree length for \acronym{RHI} calculation.
#' \insertCite{Steell2023;textual}{TreeSearch} recommend 1000, but suggest that
#' 100 may suffice.
#' If zero (the default), the \acronym{RHI} is not calculated.
#' @inheritParams CharacterLength
#' 
#' @return `Consistency()` returns a matrix with named columns specifying the 
#' consistency index (`ci`),
#' retention index (`ri`),
#' rescaled consistency index (`rc`) and
#' relative homoplasy index (`rhi`).
#' 
#' @examples 
#' data(inapplicable.datasets)
#' dataset <- inapplicable.phyData[[4]]
#' head(Consistency(dataset, TreeTools::NJTree(dataset), nRelabel = 10))
#' @references \insertAllCited{}
#' @template MRS
#' @export
Consistency <- function (dataset, tree, nRelabel = 0, compress = FALSE) {
  dsTips <- TipLabels(dataset)
  trTips <- TipLabels(tree)
  if (!setequal(dsTips, trTips)) {
    dsHas <- setdiff(dsTips, trTips)
    trHas <- setdiff(trTips, dsTips)
    stop("Tip label mismatch: ",
         if (length(dsHas)) "\n   `dataset` has ", paste(dsHas, collapse = ", "),
         if (length(trHas)) "\n   `tree` has ", paste(trHas, collapse = ", ")
    )
  }
  minLength <- MinimumLength(dataset, compress = TRUE) # Farris's m
  maxLength <- MaximumLength(dataset, compress = TRUE) # Farris's g
  tree <- Postorder(tree)
  obsLength <- CharacterLength(tree, dataset, compress = TRUE) # farris's s
  
  extra <- obsLength - minLength # Farris's h
  maxHomoplasy <- (maxLength - minLength) # g - m
  
  ci <- minLength / obsLength # Farris's c = m / s
  distortion <- extra / maxHomoplasy # Farris's d = h / (g - m)
  
  ri <- (maxLength - obsLength) / maxHomoplasy
  
  rc <- ri * minLength / obsLength

  if (nRelabel > 0) {
    medLength <- ExpectedLength(dataset, tree, nRelabel, compress = TRUE)
    expHomoplasy <- medLength - minLength
    rhi <- extra / expHomoplasy
  } else {
    rhi <- NA
  }
  
  ret <- cbind(ci = ci, ri = ri, rc = rc, rhi = rhi)
  
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(dataset, "index"), ]
  }
}


#' @importFrom fastmap fastmap
.CharLengthCache <- fastmap()

#' Expected length
#' 
#' For a given dataset and tree topology, `ExpectedLength()` estimates the
#' length expected if the states of each character are shuffled randomly
#' across the leaves.
#' 
#' @references \insertAllCited{}
#' @inheritParams Consistency
#' 
#' @return `ExpectedLength()` returns a numeric vector stating the median
#' length of each character in `dataset` on `tree` after `nRelabel` random
#' relabelling of leaves.
#' 
#' @export
#' @importFrom stats median
#' @importFrom stringi stri_paste
#' @template MRS
ExpectedLength <- function(dataset, tree, nRelabel = 1000, compress = FALSE) {
  .CheckDataCharLen(dataset)
  .CheckTreeCharLen(tree)
  tipLabel <- tree[["tip.label"]]
  tree <- .TreeForTaxa(tree, names(dataset))
  
  mat <- do.call(rbind, dataset)
  at <- attributes(dataset)
  contrast <- at[["contrast"]]
  
  rewritten <- apply(mat, 2, .SortTokens,
                     contr = apply(contrast, 1, .Bin),
                     inapp = match("-", at[["levels"]], nomatch = NA_integer_))
  rwMax <- max(rewritten)
  rwLevels <- c("-", seq_len(log2(rwMax)))
  nLevels <- length(rwLevels)
  rwContrast <- t(vapply(seq_len(rwMax), function(x) {
    as.integer(intToBits(x)[1:nLevels])
  }, integer(nLevels)))
  
  .LengthForChar <- function(x) {
    key <- stri_paste(c(nRelabel, x), collapse = ",")
    if (.CharLengthCache$has(key)) {
      .CharLengthCache$get(key)
    } else {
      patterns <- apply(unname(unique(t(
        as.data.frame(replicate(nRelabel, sample(rep(seq_along(x), x))))))),
        2, I, simplify = FALSE)
      nr <- length(patterns[[1]])
      phy <- structure(
        setNames(patterns, TipLabels(tree)),
        "weight" = rep(1, nr),
        nr = nr,
        nc = nLevels,
        index = seq_len(nr),
        levels = rwLevels,
        type = "USER",
        contrast = rwContrast,
        class = "phyDat")
      ret <- median(FastCharacterLength(tree, phy))
      .CharLengthCache$set(key, ret)
      ret
    }
  }
  
  exp <- apply(apply(rewritten, 2, tabulate, max(rewritten)), 2, .LengthForChar)
  
  # Return:
  if (compress) {
    exp
  } else {
    exp[at[["index"]]]
  }
}


.Bin <- function(x) {
  sum(2 ^ (seq_along(x)[as.logical(x)] - 1))
}


# Relabel a character such that 1 is the most common; then 2, etc.
# @param char integer vector: row of contrast matrix that applies to each taxon
# @param contr binary representation of contrast matrix
# @param inapp which level corresponds to the inapplicable state?
# @value `char` relabelled according to a binary contrast matrix in which the 
# first state corresponds to the inapplicable token
#' @importFrom fastmatch %fin% fmatch
.SortTokens <- function(char, contr, inapp = NA_integer_) {
  if (is.na(inapp)) {
    # Add a dummy inapplicable token
    inapp <- 2 ^ (floor(max(log2(contr[char]))) + 1)
    contr <- c(contr, inapp)
  }
  logCont <- log2(contr)
  maxToken <- floor(max(logCont))
  if (maxToken > 32) {
    # Too big for intToBits
    stop("This many tokens are not supported; contact the maintainer for help")
  }
  nWhole <- 2 ^ maxToken
  maxN <- nWhole + nWhole - 1 # Token that is completely ambiguous
  ambig <- logCont != floor(logCont) # Is contrast entry (partly) ambiguous?
  inappToken <- contr == inapp
  tokensToSort <- contr[!inappToken & !ambig]
  tab <- tabulate(contr[char], maxN)
  # mapping maps each token or set of tokens to its new label
  mapping <- integer(maxN)
  mapping[[inapp]] <- 1
  mapping[tokensToSort[order(tab[tokensToSort], decreasing = TRUE)]] <- 
    2 ^ seq_along(tokensToSort)
  
  nAssigned <- log2(nWhole) + 1
  wholes <- mapping[2 ^ (seq_len(nAssigned) - 1)]
  
  ambigTokens <- contr[ambig & seq_along(contr) %fin% char]
  mapping[ambigTokens] <- apply(matrix(as.logical(intToBits(contr[ambig])), 32),
                                2, function(x) sum(wholes[x]))
  
  # Return:
  mapping[contr[char]]
}
