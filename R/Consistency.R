#' Consistency and retention "indices"
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
#' \insertCite{Archie1989,Brooks1986,Steell2025}{TreeSearch}.
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
#' The **relative homoplasy index** \insertCite{Steell2025}{TreeSearch} is
#' the ratio of the observed excess tree length to the excess tree length
#' due to chance, taken as the median score of a character when the leaves
#' of the given tree are randomly shuffled.
#' 
#' The lengths of characters including inapplicable tokens are calculated
#' following \insertCite{Brazeau2019;textual}{TreeSearch}, matching their
#' default treatment in [`TreeLength()`].
#' 
#' @param nRelabel Integer specifying how many times to relabel leaves when
#' computing MCMC estimate of null tree length for \acronym{RHI} calculation.
#' \insertCite{Steell2025;textual}{TreeSearch} recommend 1000, but suggest that
#' 100 may suffice.
#' If zero (the default), the \acronym{RHI} is not calculated.
#' @inheritParams CharacterLength
#' 
#' @return `Consistency()` returns a matrix with named columns specifying the
#' consistency index (`ci`),
#' retention index (`ri`),
#' rescaled consistency index (`rc`) and
#' relative homoplasy index (`rhi`).
#' `ci` is `NaN` for a constant character, for which both the observed and
#' minimum length are zero.
#' `ri` and `rc` are `NaN` when the maximum and minimum length coincide, as
#' for a constant or an autapomorphic character.
#' `rhi` is `NA` throughout if `nRelabel = 0`, as it is then not calculated.
#' Otherwise `rhi` is `NaN` when the observed length already equals the
#' minimum length and the median length under random leaf relabelling also
#' equals the minimum; if only the median length equals the minimum, `rhi`
#' is `Inf`.
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
    ret[attr(dataset, "index"), , drop = FALSE]
  }
}


.CharLengthCache <- new.env(hash = TRUE, parent = emptyenv())

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
#' @family tree scoring
#' @template MRS
ExpectedLength <- function(dataset, tree, nRelabel = 1000, compress = FALSE) {
  .CheckDataCharLen(dataset)
  .CheckTreeCharLen(tree)
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
  
  # Key on the unlabelled rooted shape, which is what the sampled distribution
  # is a function of: leaf states are permuted uniformly, and relabelling
  # composes with a uniform permutation to leave it uniform, so any two trees
  # of the same shape are sampling the same distribution.  Keying on the
  # labelled topology instead would be sound but strictly weaker -- identical
  # topologies are a subset of identical shapes, so it would miss every reuse
  # this catches and none of its own.  Rooting is part of the shape, as these
  # characters may contain inapplicable tokens, whose lengths are not
  # rooting-invariant.
  treeKey <- .ShapeKey(tree)
  # Cache per shape, and within that per character, rather than pasting both
  # into one key: that keeps the shape key out of every character's entry, and
  # leaves no ambiguity about where the shape key ends and the counts begin.
  treeCache <- .CharLengthCache[[treeKey]]
  if (is.null(treeCache)) {
    treeCache <- new.env(hash = TRUE, parent = emptyenv())
    .CharLengthCache[[treeKey]] <- treeCache
  }

  .LengthForChar <- function(x) {
    key <- paste(c(nRelabel, x), collapse = ",")
    if (!is.null(treeCache[[key]])) {
      treeCache[[key]]
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
      treeCache[[key]] <- ret
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


# Canonical identifier of a rooted tree's unlabelled shape, after
# Aho, Hopcroft & Ullman: a leaf encodes as `01`, and an internal node wraps
# its children's codes, sorted into a fixed order, in `0`...`1`.  Sorting is
# what makes the code canonical, so it is already invariant to edge order and
# to node rotation, and two rooted shapes are isomorphic exactly if their codes
# agree.  Unlike `TreeTools::RootedTreeShape()`, which enumerates shapes into
# an integer and so stops at 55 leaves, this is bounded only by string length.
# @param tree A rooted, binary tree of class `phylo`.
# @return A string identifying the shape of `tree`.
#' @importFrom TreeTools NTip Postorder
.ShapeKey <- function(tree) {
  edge <- Postorder(tree)[["edge"]]
  nTip <- NTip(tree)
  code <- character(max(edge))
  code[seq_len(nTip)] <- "01"
  kids <- vector("list", max(edge))
  # Postorder guarantees that a node's children are coded before the edge that
  # subtends it is read, so a single pass suffices.
  for (i in seq_len(dim(edge)[[1]])) {
    parent <- edge[[i, 1]]
    kids[[parent]] <- c(kids[[parent]], code[[edge[[i, 2]]]])
    if (length(kids[[parent]]) == 2L) {
      code[[parent]] <- paste0("0", paste(sort(kids[[parent]],
                                               method = "radix"),
                                          collapse = ""), "1")
    }
  }
  bits <- as.integer(strsplit(code[[edge[[dim(edge)[[1]], 1]]]], "",
                              fixed = TRUE)[[1]]) == 1L
  # Pack to bytes for compactness.  Padding to a byte boundary could otherwise
  # conflate shapes whose codes differ only in length, so the leaf count leads.
  bits <- c(bits, rep(FALSE, (-length(bits)) %% 8))
  paste0(nTip, ":", paste(as.character(packBits(bits, "raw")), collapse = ""))
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
  wholeBits <- 2 ^ (seq_len(nAssigned) - 1)
  # A state that never occurs on its own -- only ever within an ambiguous
  # (polymorphic) token -- has no row in `mapping` yet.  Give it its own
  # unused code so that ambiguous tokens referencing it still sum to a
  # meaningful value, rather than silently contributing zero.
  unassigned <- wholeBits[mapping[wholeBits] == 0]
  if (length(unassigned)) {
    mapping[unassigned] <- 2 ^ (length(tokensToSort) + seq_along(unassigned))
  }
  wholes <- mapping[wholeBits]

  ambigTokens <- contr[ambig & seq_along(contr) %fin% char]
  mapping[ambigTokens] <- apply(matrix(as.logical(intToBits(ambigTokens)), 32),
                                2, function(x) sum(wholes[x]))
  
  # Return:
  mapping[contr[char]]
}
