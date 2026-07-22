#' Define character hierarchy for inapplicable data
#'
#' Specify the dependency structure between characters in a morphological
#' dataset that uses reductive coding.  A "controlling primary" character
#' (typically presence/absence of a structure) determines whether its
#' associated "secondary" characters are applicable.  Secondary characters
#' can in turn control tertiary characters, and so on.
#'
#' This hierarchy is required by the HSJ
#' \insertCite{Hopkins2021}{TreeSearch} and step-matrix
#' \insertCite{Goloboff2021}{TreeSearch} approaches to inapplicable
#' characters, and is passed to [`MaximizeParsimony()`] via the `hierarchy`
#' argument.
#'
#' @param ... Named arguments where each name is the index of a controlling
#'   character (coerced to integer) and each value is an integer vector of
#'   the character indices it controls.  Use nested [`list()`]s for deeper
#'   hierarchies (see Examples).
#'
#' @return An object of class `"CharacterHierarchy"`.
#'
#' @examples
#' # Simple: character 1 controls characters 2-5
#' h <- CharacterHierarchy("1" = 2:5)
#'
#' # Multiple controlling primaries
#' h <- CharacterHierarchy("1" = 2:5, "6" = 7:8)
#'
#' # Nested: char 1 controls 2-5; char 3 further controls 9-10
#' h <- CharacterHierarchy("1" = list(2, 3, 4, 5, "3" = 9:10))
#'
#' @references
#' \insertAllCited{}
#' @family tree scoring
#' @seealso [MaximizeParsimony()], [HierarchyFromNames()]
#' @export
CharacterHierarchy <- function(...) {
  args <- list(...)
  if (length(args) == 0L) {
    stop("At least one controlling character must be specified.")
  }
  tree <- .ParseHierarchyArgs(args)
  structure(tree, class = "CharacterHierarchy")
}

# Parse user args into a normalized tree structure.
# Returns a list of nodes, each:
#   list(controlling = int, dependents = int[], children = list(<node>, ...))
# "children" are sub-hierarchies (controlling secondaries).
.ParseHierarchyArgs <- function(args) {
  if (is.null(names(args)) || any(names(args) == "")) {
    stop("Every element of `...` must be named with the controlling ",
         "character index.")
  }
  controllingIndices <- suppressWarnings(as.integer(names(args)))
  if (anyNA(controllingIndices)) {
    stop("Controlling character names must be integer indices.")
  }

  lapply(seq_along(args), function(i) {
    ctrl <- controllingIndices[i]
    val <- args[[i]]
    .ParseOneBlock(ctrl, val)
  })
}

# Parse a single controlling-character block.
# val can be:
#   - integer vector: simple list of dependent character indices
#   - list with mixed named/unnamed elements: unnamed = dependents,
#     named = sub-hierarchies (controlling secondaries)
.ParseOneBlock <- function(ctrl, val) {
  if (is.numeric(val) && is.null(names(val))) {
    # Simple case: vector of dependent indices
    return(list(
      controlling = as.integer(ctrl),
      dependents = as.integer(val),
      children = list()
    ))
  }
  if (is.list(val)) {
    nms <- names(val)
    if (is.null(nms)) nms <- rep("", length(val))
    dependents <- integer(0)
    children <- list()
    for (j in seq_along(val)) {
      if (nms[j] == "") {
        # Unnamed: a dependent character index
        dependents <- c(dependents, as.integer(val[[j]]))
      } else {
        # Named: a sub-hierarchy
        subCtrl <- suppressWarnings(as.integer(nms[j]))
        if (is.na(subCtrl)) {
          stop("Sub-hierarchy names must be integer character indices, got '",
               nms[j], "'.")
        }
        # The sub-controlling character is also a dependent of this block
        dependents <- c(dependents, subCtrl)
        children <- c(children, list(.ParseOneBlock(subCtrl, val[[j]])))
      }
    }
    return(list(
      controlling = as.integer(ctrl),
      # A sub-controller may also be listed as an explicit dependent (e.g.
      # `list(2, 3, 4, 5, "3" = 9:10)`); keep it once so ValidateHierarchy()
      # does not flag it as appearing in multiple blocks.
      dependents = unique(dependents),
      children = children
    ))
  }
  # Scalar
  list(
    controlling = as.integer(ctrl),
    dependents = as.integer(val),
    children = list()
  )
}

#' @export
print.CharacterHierarchy <- function(x, ...) {
  cat("CharacterHierarchy\n")
  .PrintBlock <- function(node, indent = 1L) {
    pad <- strrep("  ", indent)
    leafDeps <- setdiff(
      node$dependents,
      vapply(node$children, `[[`, integer(1), "controlling")
    )
    cat(sprintf("%sChar %d controls: {%s}\n",
                pad, node$controlling,
                paste(node$dependents, collapse = ", ")))
    for (child in node$children) {
      .PrintBlock(child, indent + 1L)
    }
  }
  for (node in x) {
    .PrintBlock(node)
  }
  invisible(x)
}

#' Validate a CharacterHierarchy against a dataset
#'
#' Check that a [`CharacterHierarchy`] object is consistent with a
#' [`phyDat`][phangorn::phyDat] dataset: character indices exist,
#' controlling characters are binary (absent/present), secondaries are
#' coded inapplicable where expected, and no character appears in
#' multiple blocks.
#'
#' @param hierarchy A [`CharacterHierarchy`] object.
#' @param dataset A `phyDat` object.
#'
#' @return `hierarchy`, invisibly (called for side effects: stops with an
#'   informative error if validation fails).
#'
#' @keywords internal
#' @importFrom utils head
#' @export
ValidateHierarchy <- function(hierarchy, dataset) {
  if (!inherits(hierarchy, "CharacterHierarchy")) {
    stop("`hierarchy` must be a CharacterHierarchy object.")
  }
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }

  nChar <- length(attr(dataset, "index"))
  allLevels <- attr(dataset, "allLevels")
  levels <- attr(dataset, "levels")
  contrast <- attr(dataset, "contrast")

  # Identify the inapplicable token
  inappToken <- "-"
  if (!inappToken %in% allLevels) {
    stop("Dataset does not contain an inapplicable token ('-').")
  }

  # Build the original character matrix
  idx <- attr(dataset, "index")
  origMat <- do.call(rbind, lapply(dataset, function(x) {
    allLevels[x[idx]]
  }))

  # Identify the "0" state (absence) in the controlling primary
  absenceState <- "0"

  # Track all characters claimed by any block

  claimed <- integer(0)

  .ValidateBlock <- function(node, depth = 1L) {
    ctrl <- node$controlling
    deps <- node$dependents

    # Check indices exist
    allIdx <- c(ctrl, deps)
    bad <- allIdx[allIdx < 1L | allIdx > nChar]
    if (length(bad) > 0L) {
      stop(sprintf(
        "Character index(es) %s out of range [1, %d].",
        paste(bad, collapse = ", "), nChar
      ))
    }

    # Check no double-claiming
    overlap <- intersect(allIdx, claimed)
    if (length(overlap) > 0L) {
      stop(sprintf(
        "Character(s) %s appear in multiple hierarchy blocks.",
        paste(overlap, collapse = ", ")
      ))
    }
    claimed <<- c(claimed, allIdx)

    # Check controlling character is binary (has exactly states "0" and "1",
    # possibly with inapplicable/missing)
    ctrlVals <- unique(origMat[, ctrl])
    ctrlInformative <- setdiff(ctrlVals, c("?", "-"))
    if (!all(ctrlInformative %in% c("0", "1"))) {
      stop(sprintf(
        paste0("Controlling character %d must be binary (states '0' and '1'),",
               " but has states: %s."),
        ctrl, paste(ctrlInformative, collapse = ", ")
      ))
    }

    # Check secondaries are "-" where controlling is "0"
    absentTaxa <- which(origMat[, ctrl] == absenceState)
    if (length(absentTaxa) > 0L) {
      for (d in deps) {
        depVals <- origMat[absentTaxa, d]
        badTaxa <- which(!depVals %in% c("-", "?"))
        if (length(badTaxa) > 0L) {
          badNames <- rownames(origMat)[absentTaxa[badTaxa]]
          stop(sprintf(
            paste0("Secondary character %d has non-inapplicable values for ",
                   "taxa where controlling character %d is absent: %s."),
            d, ctrl, paste(head(badNames, 5), collapse = ", ")
          ))
        }
      }
    }

    # Recurse into children
    for (child in node$children) {
      .ValidateBlock(child, depth + 1L)
    }
  }

  for (node in hierarchy) {
    .ValidateBlock(node)
  }

  invisible(hierarchy)
}


#' Construct a CharacterHierarchy from TNT-style character names
#'
#' Parse character names following the TNT convention where controlling
#' characters are named `sup_<tag>` and their dependent characters are
#' named `sub_<tag>[_suffix]`.  Tags must match between a controlling
#' character and its dependents.  Nested hierarchies are detected when a
#' `sub_` character is also a `sup_` for further characters.
#'
#' @param charNames Character vector of names, one per original character.
#'
#' @return A [`CharacterHierarchy`] object, or `NULL` if no hierarchy is
#'   detected.
#'
#' @examples
#' names <- c("sup_tail", "sub_tail_colour", "sub_tail_shape",
#'             "sup_wing", "sub_wing_venation", "eyes")
#' HierarchyFromNames(names)
#'
#' @family tree scoring
#' @seealso [CharacterHierarchy()]
#' @export
HierarchyFromNames <- function(charNames) {
  if (!is.character(charNames) || length(charNames) == 0L) {
    stop("`charNames` must be a non-empty character vector.")
  }

  # Find sup_ and sub_ characters
  supIdx <- grep("^sup_", charNames)
  subIdx <- grep("^sub_", charNames)

  if (length(supIdx) == 0L) {
    return(NULL)
  }

  # Extract tags
  supTags <- sub("^sup_", "", charNames[supIdx])
  subTagsFull <- sub("^sub_", "", charNames[subIdx])
  # The tag is the first component before any additional underscore-suffix
  # e.g. "sub_tail_colour" → tag = "tail"
  subTags <- sub("_.*", "", subTagsFull)

  # Build mapping: tag → controlling index, tag → dependent indices
  tagToSup <- setNames(supIdx, supTags)

  # Group sub characters by tag
  tagToSubs <- split(subIdx, subTags)

  # Check for sub_ characters referencing nonexistent sup_ tags
  orphanTags <- setdiff(names(tagToSubs), supTags)
  if (length(orphanTags) > 0L) {
    warning(sprintf(
      "sub_ characters reference tags with no corresponding sup_: %s",
      paste(orphanTags, collapse = ", ")
    ))
  }

  # Detect nested hierarchies: a sub_ character that is also a sup_
  # Find sub_ chars that are also in supIdx
  subAlsoSup <- intersect(subIdx, supIdx)

  # Build hierarchy
  # First pass: create flat blocks for all sup_ tags
  args <- list()
  for (tag in supTags) {
    ctrl <- tagToSup[[tag]]
    subs <- tagToSubs[[tag]]
    if (is.null(subs)) subs <- integer(0)

    # Check which subs are themselves controlling (nested hierarchy)
    nestedSubs <- intersect(subs, supIdx)
    flatSubs <- setdiff(subs, supIdx)

    if (length(nestedSubs) == 0L) {
      # Simple block
      args[[as.character(ctrl)]] <- as.integer(subs)
    } else {
      # Nested: build list with named sub-hierarchies
      block <- as.list(as.integer(flatSubs))
      for (ns in nestedSubs) {
        nsTag <- supTags[supIdx == ns]
        nsSubs <- tagToSubs[[nsTag]]
        if (is.null(nsSubs)) nsSubs <- integer(0)
        block[[as.character(ns)]] <- as.integer(nsSubs)
      }
      args[[as.character(ctrl)]] <- block
    }
  }

  # Filter out sup_ chars whose index also appears in subIdx
  # (they'll be included as children of their parent)
  topLevelSup <- setdiff(supIdx, subIdx)
  if (length(topLevelSup) == 0L) {
    # All sup_ characters are also sub_ — circular or all nested.
    # Fall back to treating all as top-level with a warning.
    warning("All sup_ characters are also sub_ characters. ",
            "Treating all as top-level.")
    topLevelSup <- supIdx
  }
  topLevelCtrls <- as.character(topLevelSup)
  args <- args[topLevelCtrls]

  do.call(CharacterHierarchy, args)
}


#' Extract all character indices from a hierarchy
#'
#' Returns all character indices (controlling + dependent) referenced by
#' a [`CharacterHierarchy`], useful for partitioning characters into
#' hierarchy vs. non-hierarchy sets.
#'
#' @param hierarchy A [`CharacterHierarchy`] object.
#'
#' @return An integer vector of character indices (unsorted, may contain
#'   duplicates if the hierarchy is malformed).
#'
#' @keywords internal
#' @export
HierarchyChars <- function(hierarchy) {
  .CollectIndices <- function(node) {
    c(node$controlling, node$dependents,
      unlist(lapply(node$children, .CollectIndices)))
  }
  unique(unlist(lapply(hierarchy, .CollectIndices)))
}


#' List top-level controlling characters
#'
#' @param hierarchy A [`CharacterHierarchy`] object.
#' @return Integer vector of top-level controlling character indices.
#' @keywords internal
#' @export
HierarchyControlling <- function(hierarchy) {
  vapply(hierarchy, `[[`, integer(1), "controlling")
}


# Build the tip-labels matrix for HSJ scoring.
#
# Converts a phyDat dataset into an integer matrix of per-tip, per-character
# state labels (0-based) for the C++ HSJ scorer: length(dataset) rows (tips) by
# length(attr(dataset, "index")) columns (original characters).
.BuildTipLabels <- function(dataset) {
  idx <- attr(dataset, "index")
  nTip <- length(dataset)
  nChar <- length(idx)

  # dataset is a list of integer vectors (pattern indices per tip)
  # Expand via index to original characters, convert to 0-based
  mat <- matrix(0L, nrow = nTip, ncol = nChar)
  for (t in seq_len(nTip)) {
    patternTokens <- dataset[[t]]    # token indices for each pattern
    mat[t, ] <- patternTokens[idx] - 1L  # 0-based
  }
  mat
}


# Identify the primary "absent" state for HSJ scoring.
#
# Returns the 0-based token index of the controlling primary character's
# *absent* state, for the C++ HSJ scorer's `absent_state` argument.
#
# Under reductive coding (Hopkins & St John 2021) the primary codes a
# structure's presence/absence, conventionally "0" = absent, "1" = present.
# The index of "0" depends on the dataset's `levels` ordering (e.g. it is 1 for
# c("-", "0", "1") but 0 for c("0", "1")), so it must be computed rather than
# hard-coded.  The inapplicable token "-" is also treated as absent by the
# scorer; if no "0" level exists, the index of "-" is returned.
.HSJAbsentState <- function(dataset) {
  lv <- attr(dataset, "levels")
  idx <- match("0", lv)
  if (is.na(idx)) {
    idx <- match("-", lv)
  }
  if (is.na(idx)) 0L else as.integer(idx - 1L)
}


# Convert a CharacterHierarchy into a flat list of hierarchy blocks for the C++
# ts_hsj_score() bridge.  Each block is a list with `primary` (0-based) and
# `secondaries` (0-based integer vector).
.HierarchyToBlocks <- function(hierarchy) {
  .FlattenBlock <- function(node) {
    block <- list(
      primary = node$controlling - 1L,
      secondaries = node$dependents - 1L
    )
    childBlocks <- lapply(node$children, .FlattenBlock)
    c(list(block), unlist(childBlocks, recursive = FALSE))
  }
  unlist(lapply(hierarchy, .FlattenBlock), recursive = FALSE)
}


# Compute non-hierarchy pattern weights: given a phyDat dataset and a
# CharacterHierarchy, return the integer weight vector (same length as
# attr(dataset, "weight")) with hierarchy characters' contributions subtracted.
# Patterns appearing only in hierarchy characters end up with weight 0.
.NonHierarchyWeights <- function(dataset, hierarchy) {
  w <- attr(dataset, "weight")
  idx <- attr(dataset, "index")
  hChars <- HierarchyChars(hierarchy)

  adjusted <- as.integer(w)
  for (ci in hChars) {
    if (ci < 1L || ci > length(idx)) next
    pat <- idx[ci]
    if (pat >= 1L && pat <= length(adjusted) && adjusted[pat] > 0L) {
      adjusted[pat] <- adjusted[pat] - 1L
    }
  }
  adjusted
}


# Generate resampled weights for hierarchical resampling.
#
# Instead of treating every character independently, groups characters into
# resampling units: each non-hierarchy character is one unit, and each
# top-level hierarchy block (primary + all dependents, recursively) is one
# unit.  Jackknife or bootstrap operates on these units.
#
# Returns a list with:
#   nonHierarchyWeights: pattern weights for Fitch scoring (non-hierarchy
#     chars only, reflecting which free chars were sampled)
#   blockCounts: integer vector (length = number of top-level blocks)
#     giving how many times each block was sampled (0/1 for jackknife,
#     0+ for bootstrap)
.HierarchicalResampleWeights <- function(dataset, hierarchy, bootstrap,
                                         proportion) {
  idx <- attr(dataset, "index")
  nPatterns <- length(attr(dataset, "weight"))
  nChars <- length(idx)

  # Collect chars per top-level block (includes nested dependents)
  .CollectAll <- function(node) {
    c(node$controlling, node$dependents,
      unlist(lapply(node$children, .CollectAll)))
  }
  nBlocks <- length(hierarchy)
  blockChars <- lapply(hierarchy, function(node) unique(.CollectAll(node)))
  hCharsSet <- unique(unlist(blockChars))

  freeChars <- setdiff(seq_len(nChars), hCharsSet)
  nFree <- length(freeChars)
  nUnits <- nFree + nBlocks

  if (nUnits < 2L) {
    # Degenerate: can't jackknife with < 2 units
    return(list(
      nonHierarchyWeights = .NonHierarchyWeights(dataset, hierarchy),
      blockCounts = rep(1L, nBlocks)
    ))
  }

  if (bootstrap) {
    sampled <- sample.int(nUnits, nUnits, replace = TRUE)
  } else {
    nKeep <- max(1L, ceiling(proportion * nUnits))
    nKeep <- min(nKeep, nUnits - 1L)
    sampled <- sample.int(nUnits, nKeep, replace = FALSE)
  }

  unitCounts <- tabulate(sampled, nbins = nUnits)

  # Non-hierarchy pattern weights from retained free chars
  nhWeights <- integer(nPatterns)
  for (i in seq_len(nFree)) {
    if (unitCounts[i] > 0L) {
      pat <- idx[freeChars[i]]
      nhWeights[pat] <- nhWeights[pat] + unitCounts[i]
    }
  }

  blockCounts <- unitCounts[nFree + seq_len(nBlocks)]

  list(
    nonHierarchyWeights = nhWeights,
    blockCounts = blockCounts
  )
}
