#' Calculate the parsimony score of a tree given a dataset
#'
#' `TreeLength()` calculates a parsimony score for a tree.
#' Trees may be scored using equal weights, implied weights
#' \insertCite{Goloboff1993}{TreeSearch}, or profile parsimony
#' \insertCite{Faith2001}{TreeSearch}.
#' Inapplicable characters are handled using the algorithm of
#' \insertCite{Brazeau2019;textual}{TreeSearch} by default, or
#' alternatively using the hierarchical scoring of
#' \insertCite{Hopkins2021;textual}{TreeSearch} when
#' `inapplicable = "hsj"` and a [`CharacterHierarchy`] is provided.
#'
#' @param tree A tree of class `phylo`, a list thereof (optionally of class
#' `multiPhylo`), or an integer -- in which case `tree` random trees will be 
#' uniformly sampled.
#' @inheritParams MaximizeParsimony
#' 
#' @return `TreeLength()` returns a numeric vector containing the score for
#' each tree in `tree`.
#' 
#' @examples
#' data("inapplicable.datasets")
#' tree <- TreeTools::BalancedTree(inapplicable.phyData[[1]])
#' TreeLength(tree, inapplicable.phyData[[1]])
#' TreeLength(tree, inapplicable.phyData[[1]], concavity = 10)
#' \donttest{ # PrepareDataProfile() and random-tree scoring are slower:
#' TreeLength(tree, inapplicable.phyData[[1]], concavity = "profile")
#' TreeLength(5, inapplicable.phyData[[1]])
#'
#' # HSJ scoring with a character hierarchy
#' dataset6 <- inapplicable.phyData[["Vinther2008"]]
#' hier <- CharacterHierarchy("1" = 2:3)
#' tree6 <- TreeTools::BalancedTree(dataset6)
#' TreeLength(tree6, dataset6, hierarchy = hier, inapplicable = "hsj")
#' }
#' @seealso 
#' - Conduct tree search using [`MaximizeParsimony()`] (command line) or
#' [`EasyTrees()`] (graphical user interface).
#' 
#' - See score for each character: [`CharacterLength()`].
#' @family tree scoring 
#' 
#' @references
#' \insertAllCited{}
#' @author Martin R. Smith
#' @importFrom fastmatch %fin%
#' @importFrom TreeTools Renumber RenumberTips TreeIsRooted
#' @export
TreeLength <- function(tree, dataset, concavity = Inf,
                       extended_iw = TRUE,
                       xpiwe_r = 0.5,
                       xpiwe_max_f = 5,
                       hierarchy = NULL, inapplicable = "bgs",
                       hsj_alpha = 1.0) {
  UseMethod("TreeLength")
}

#' @rdname TreeLength
#' @export
TreeLength.phylo <- function(tree, dataset, concavity = Inf,
                              extended_iw = TRUE,
                              xpiwe_r = 0.5,
                              xpiwe_max_f = 5,
                              hierarchy = NULL, inapplicable = "bgs",
                              hsj_alpha = 1.0) {
  tipLabels <- tree[["tip.label"]]
  
  if (!TreeIsRooted(tree)) {
    stop("`tree` must be rooted; try RootTree(tree)")
  }
  
  nTip <- length(tipLabels)
  edge <- tree[["edge"]]
  if (dim(edge)[1] != nTip + nTip - 2) {
    stop("`tree` must be binary")
  }
  
  if (length(setdiff(tipLabels, names(dataset)))) {
      stop("Missing in `dataset`: ",
           paste(setdiff(tipLabels, names(dataset)), collapse = ", "))
  }
  
  if (is.null(attr(dataset, "levels")) || ncol(attr(dataset, "contrast")) == 0L) {
    return(0L)
  }

  if (nTip < length(dataset)) {
    dataset <- .Recompress(dataset[tree[["tip.label"]]])
  }

  # --- Validate inapplicable-handling parameters ---
  inapplicable <- tolower(inapplicable)
  if (inapplicable == "brazeau") inapplicable <- "bgs"
  inapplicable <- match.arg(inapplicable, c("bgs", "hsj", "xform", "missing"))
  if (inapplicable == "missing") {
    # Gaps as missing data: recode, then score with the standard Fitch engine.
    # Placed after .Recompress() (whose round-trip would otherwise restore the
    # "-" level); PrepareDataIW() preserves the recode and profile parsimony
    # handles gaps itself.
    dataset <- .GapsAsMissing(dataset)
    inapplicable <- "bgs"
  }
  useHSJ <- !is.null(hierarchy) && identical(inapplicable, "hsj")
  if (inapplicable != "bgs") {
    if (is.null(hierarchy)) {
      stop("A `hierarchy` is required when inapplicable = \"", inapplicable,
           "\". See ?CharacterHierarchy.")
    }
    if (!inherits(hierarchy, "CharacterHierarchy")) {
      stop("`hierarchy` must be a CharacterHierarchy object.")
    }
    ValidateHierarchy(hierarchy, dataset)
    if (.UseProfile(concavity)) {
      stop("Profile parsimony is not currently supported with inapplicable = \"",
           inapplicable, "\".")
    }
    if (is.finite(concavity)) {
      stop("Implied weighting is not currently supported with inapplicable = \"",
           inapplicable, "\".")
    }
  }
  useXform <- !is.null(hierarchy) && identical(inapplicable, "xform")
  if (!is.numeric(hsj_alpha) || length(hsj_alpha) != 1L ||
      hsj_alpha < 0 || hsj_alpha > 1) {
    stop("`hsj_alpha` must be a single number in [0, 1].")
  }

  if (is.finite(concavity)) {
    if (concavity <= 0) {
      stop("`concavity` must be positive (or Inf for equal weights, ",
           "or \"profile\" for profile parsimony).")
    }
    if (!("min.length" %fin% names(attributes(dataset)))) {
      dataset <- PrepareDataIW(dataset)
    }
    at <- attributes(dataset)
    nChar  <- at[["nr"]] # strictly, transformation series patterns
                    # these'll be upweighted later
    weight <- at[["weight"]]
    steps <- CharacterLength(tree, dataset, compress = TRUE)
    minLength <- at[["min.length"]]
    homoplasies <- steps - minLength
    
    # This check was once triggered - possibly fixed but remains
    # under investigation...
    if (any(homoplasies < 0)) { #nocov start
      stop("Minimum steps have been miscalculated.\n", 
           "       Please report this bug at:\n", 
           "       https://github.com/ms609/TreeSearch/issues/new\n\n",
           "       See above for full tree: ", dput(tree))
    } #nocov end
    if (isTRUE(extended_iw)) {
      obsCount <- .ObsCount(dataset)
      nTaxa <- length(dataset)
      # Goloboff (2014) Extension 3, verified against TNT 1.6:
      # f = 1 + r * missing / obs  (NOT r * total / obs)
      f <- pmin(pmax(1 + xpiwe_r * (nTaxa - obsCount) / obsCount, 1),
                xpiwe_max_f)
      eff_k <- concavity / f
      phi <- (1 + eff_k) / (1 + concavity)
    } else {
      eff_k <- concavity
      phi <- 1
    }
    fit <- homoplasies / (homoplasies + eff_k)
    # Return:
    sum(fit * weight * phi)
    
  } else if (.UseProfile(concavity)) {
    dataset <- PrepareDataProfile(dataset)
    steps <- CharacterLength(tree, dataset, compress = TRUE)
    info <- attr(dataset, "info.amounts")
    
    # Return:
    sum(vapply(which(steps > 0), function(i) info[steps[i], i],
               double(1)) * attr(dataset, "weight")[steps > 0])
  } else if (useHSJ) {
    tree <- RenumberTips(Renumber(tree), names(dataset))
    at <- attributes(dataset)
    contrast <- at$contrast
    tip_data <- matrix(unlist(dataset, use.names = FALSE),
                       nrow = length(dataset), byrow = TRUE)
    adj_weight <- .NonHierarchyWeights(dataset, hierarchy)
    ts_hsj_score(tree[["edge"]], contrast, tip_data,
                 as.integer(adj_weight), at$levels,
                 .HierarchyToBlocks(hierarchy),
                 as.double(hsj_alpha),
                 .BuildTipLabels(dataset),
                 .HSJAbsentState(dataset))
  } else if (useXform) {
    tree <- RenumberTips(Renumber(tree), names(dataset))
    at <- attributes(dataset)
    contrast <- at$contrast
    tip_data <- matrix(unlist(dataset, use.names = FALSE),
                       nrow = length(dataset), byrow = TRUE)
    adj_weight <- as.integer(.NonHierarchyWeights(dataset, hierarchy))
    recoded <- RecodeHierarchy(dataset, hierarchy)
    xform <- .PrepareXformArgs(recoded, length(dataset))
    fitch_part <- ts_fitch_score(tree[["edge"]], contrast, tip_data,
                                 adj_weight, at$levels)
    res <- ts_sankoff_test(tree[["edge"]], xform$n_states,
                           xform$cost_matrices, xform$tip_states,
                           xform$forced_root)
    fitch_part + res$score
  } else {
    tree <- RenumberTips(Renumber(tree), names(dataset))
    at <- attributes(dataset)
    contrast <- at$contrast
    tip_data <- matrix(unlist(dataset, use.names = FALSE),
                       nrow = length(dataset), byrow = TRUE)
    ts_fitch_score(tree[["edge"]], contrast, tip_data,
                   .ScaleWeight(at$weight), at$levels)
  }
}


#' @rdname TreeLength
#' @importFrom TreeTools RandomTree
#' @export
TreeLength.numeric <- function(tree, dataset, concavity = Inf,
                               extended_iw = TRUE,
                               xpiwe_r = 0.5,
                               xpiwe_max_f = 5,
                               hierarchy = NULL, inapplicable = "bgs",
                               hsj_alpha = 1.0) {
  TreeLength(lapply(!logical(tree), RandomTree, tips = dataset), 
             dataset = dataset, concavity = concavity,
             extended_iw = extended_iw,
             xpiwe_r = xpiwe_r, xpiwe_max_f = xpiwe_max_f,
             hierarchy = hierarchy, inapplicable = inapplicable,
             hsj_alpha = hsj_alpha)
}

#' @rdname TreeLength
#' @export
TreeLength.list <- function(tree, dataset, concavity = Inf,
                            extended_iw = TRUE,
                            xpiwe_r = 0.5,
                            xpiwe_max_f = 5,
                            hierarchy = NULL, inapplicable = "bgs",
                            hsj_alpha = 1.0) {
  iw <- is.finite(concavity)
  useProfile <- .UseProfile(concavity)

  # --- Validate inapplicable-handling parameters ---
  inapplicable <- tolower(inapplicable)
  if (inapplicable == "brazeau") inapplicable <- "bgs"
  inapplicable <- match.arg(inapplicable, c("bgs", "hsj", "xform", "missing"))
  # Gaps as missing data: the recode is deferred to after any .Recompress()
  # round-trip below (which would otherwise restore the "-" level).
  useMissing <- inapplicable == "missing"
  if (useMissing) inapplicable <- "bgs"
  useHSJ <- !is.null(hierarchy) && identical(inapplicable, "hsj")
  if (inapplicable != "bgs") {
    if (is.null(hierarchy)) {
      stop("A `hierarchy` is required when inapplicable = \"", inapplicable,
           "\". See ?CharacterHierarchy.")
    }
    if (!inherits(hierarchy, "CharacterHierarchy")) {
      stop("`hierarchy` must be a CharacterHierarchy object.")
    }
    ValidateHierarchy(hierarchy, dataset)
    if (useProfile) {
      stop("Profile parsimony is not currently supported with inapplicable = \"",
           inapplicable, "\".")
    }
    if (iw) {
      stop("Implied weighting is not currently supported with inapplicable = \"",
           inapplicable, "\".")
    }
  }
  useXform <- !is.null(hierarchy) && identical(inapplicable, "xform")
  if (!is.numeric(hsj_alpha) || length(hsj_alpha) != 1L ||
      hsj_alpha < 0 || hsj_alpha > 1) {
    stop("`hsj_alpha` must be a single number in [0, 1].")
  }
  if (iw && concavity <= 0) {
    stop("`concavity` must be positive (or Inf for equal weights, ",
         "or \"profile\" for profile parsimony).")
  }

  nTip <- NTip(tree)
  if (length(unique(nTip)) > 1L) {
    stop("All trees must bear the same leaves.")
  }
  nTip <- nTip[1]
  if (nTip < length(dataset)) {
    dataset <- .Recompress(dataset[TipLabels(tree[[1]])])
  }
  if (useMissing) {
    dataset <- .GapsAsMissing(dataset)
  }

  tree[] <- RenumberTips(tree, dataset)
  needRoot <- !vapply(tree, TreeIsRooted, logical(1L))
  if (any(needRoot)) warning("Unrooted tree rooted on tip 1.")
  tree[] <- lapply(tree, function(tr) if (TreeIsRooted(tr)) tr else RootTree(tr, 1))

  nEdge <- unique(vapply(tree, function(tr) dim(tr[["edge"]])[1], integer(1)))
  if (length(nEdge) > 1L) {
    stop("Trees have different numbers of edges (",
           paste0(nEdge, collapse = ", "),
           "); try collapsing polytomies?)")
  }

  if (is.null(attr(dataset, "levels")) || ncol(attr(dataset, "contrast")) == 0L) {
    return(rep(0L, length(tree)))
  }

  # Prepare dataset for C++ engine
  if (useProfile) {
    dataset <- PrepareDataProfile(dataset)
  } else if (iw) {
    if (!("min.length" %fin% names(attributes(dataset)))) {
      dataset <- PrepareDataIW(dataset)
    }
  }

  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  weight <- .ScaleWeight(at$weight)
  levels <- at$levels

  min_steps <- if (iw) as.integer(at[["min.length"]]) else integer(0)
  concavity_val <- if (iw) concavity else Inf
  infoAmounts <- if (useProfile) at$info.amounts else NULL

  # XPIWE: per-pattern observed-taxa counts
  useXpiwe <- isTRUE(extended_iw) && iw && !useProfile
  obsCount <- if (useXpiwe) .ObsCount(dataset) else integer(0)

  if (useHSJ) {
    adj_weight <- as.integer(.NonHierarchyWeights(dataset, hierarchy))
    blocks <- .HierarchyToBlocks(hierarchy)
    alpha <- as.double(hsj_alpha)
    tip_labels <- .BuildTipLabels(dataset)
    absent_state <- .HSJAbsentState(dataset)
    vapply(tree, function(tr) {
      ts_hsj_score(tr[["edge"]], contrast, tip_data, adj_weight, levels,
                   blocks, alpha, tip_labels, absent_state)
    }, double(1))
  } else if (useXform) {
    adj_weight <- as.integer(.NonHierarchyWeights(dataset, hierarchy))
    recoded <- RecodeHierarchy(dataset, hierarchy)
    xform <- .PrepareXformArgs(recoded, length(dataset))
    vapply(tree, function(tr) {
      fitch_part <- ts_fitch_score(tr[["edge"]], contrast, tip_data,
                                   adj_weight, levels)
      res <- ts_sankoff_test(tr[["edge"]], xform$n_states,
                             xform$cost_matrices, xform$tip_states,
                             xform$forced_root)
      fitch_part + res$score
    }, double(1))
  } else {
    vapply(tree, function(tr) {
      ts_fitch_score(tr[["edge"]], contrast, tip_data, weight, levels,
                     min_steps = min_steps, concavity = concavity_val,
                     infoAmounts = infoAmounts,
                     xpiwe = useXpiwe,
                     xpiwe_r = as.double(xpiwe_r),
                     xpiwe_max_f = as.double(xpiwe_max_f),
                     obs_count = obsCount)
    }, double(1))
  }
}


#' @rdname TreeLength
#' @export
TreeLength.multiPhylo <- TreeLength.list

#' @export
TreeLength.NULL <- function(tree, dataset, concavity = Inf,
                            extended_iw = TRUE,
                            xpiwe_r = 0.5,
                            xpiwe_max_f = 5,
                            hierarchy = NULL, inapplicable = "bgs",
                            hsj_alpha = 1.0) NULL

# Pack RecodeHierarchy() output into the format ts_sankoff_test() expects.
.PrepareXformArgs <- function(recoded, n_tip) {
  chars <- recoded$sankoff_chars
  n_chars <- length(chars)
  n_states <- as.integer(vapply(chars, function(ch) ch$n_states, numeric(1)))
  forced_root <- as.integer(vapply(chars, function(ch) ch$forced_root_state, numeric(1)))
  cost_matrices <- lapply(chars, function(ch) ch$cost_matrix)
  tip_states <- matrix(0L, nrow = n_tip, ncol = n_chars)
  for (i in seq_len(n_chars)) {
    tip_states[, i] <- chars[[i]]$tip_states
  }
  list(n_states = n_states, cost_matrices = cost_matrices,
       tip_states = tip_states, forced_root = forced_root)
}

#' @rdname TreeLength
#' @export
Fitch <- function(tree, dataset) {
  .Deprecated("TreeLength")
  TreeLength(tree, dataset, Inf)
}


.CheckDataCharLen <- function(dataset) {
  if (!inherits(dataset, "phyDat")) {
    stop("Dataset must be of class phyDat, not ", class(dataset), ".")
  }
}

.CheckTreeCharLen <- function(tree) {
  if (!inherits(tree, "phylo")) {
    stop("Tree must be of class phylo, not ", class(tree), ".")
  }
  if (is.null(TipLabels(tree))) {
    stop("Tree has no labels")
  }
  if (!TreeIsRooted(tree)) {
    stop("`tree` must be rooted; try RootTree(tree)")
  }
}

#' @importFrom cli cli_alert
.DataForTaxa <- function(dataset, tipLabel) {
  dataNames <- names(dataset)
  
  if (length(tipLabel) < length(dataNames)) {
    if (all(tipLabel %in% dataNames)) {
      cli_alert(paste0(
        paste0(setdiff(dataNames, tipLabel), collapse = ", "),
        " not in tree"))
      dataset <- dataset[intersect(dataNames, tipLabel)]
    } else {
      stop("Tree tips ", 
           paste(setdiff(tipLabel, dataNames), collapse = ", "),
           " not found in dataset.")
    }
  }
  dataset
}

#' @importFrom TreeTools TipLabels KeepTip Postorder RenumberTips
.TreeForTaxa <- function(tree, dataNames) {
  tipLabel <- TipLabels(tree)
  if (length(tipLabel) > length(dataNames)) {
    cli_alert(paste0(
      paste0(setdiff(tipLabel, dataNames), collapse = ", "),
      " not in `dataset`"))
    
    tree <- KeepTip(tree, dataNames)
  }
  # Morphy requires that the tree is in postorder
  tree <- RenumberTips(Postorder(tree), dataNames)
}

#' Character length
#' 
#' Homoplasy length of each character in a dataset on a specified tree.
#' 
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#' @param compress Logical specifying whether to retain the compression of a
#' `phyDat` object or to return a vector specifying to each individual
#' character, decompressed using the dataset's `index` attribute.
#'
#' @return `CharacterLength()` returns a vector listing the contribution of each
#' character to tree score, according to the algorithm of
#' \insertCite{Brazeau2018;textual}{TreeTools}.
#'
#' @examples
#' data("inapplicable.datasets")
#' dataset <- inapplicable.phyData[[12]]
#' tree <- TreeTools::NJTree(dataset)
#' CharacterLength(tree, dataset)
#' CharacterLength(tree, dataset, compress = TRUE)
#' @template MRS
#' @family tree scoring
#' @references
#' \insertAllCited{}
#' @export
CharacterLength <- function(tree, dataset, compress = FALSE) {
  .CheckDataCharLen(dataset)
  .CheckTreeCharLen(tree)
  tipLabel <- tree[["tip.label"]]
  dataset <- .DataForTaxa(dataset, tipLabel)
  tree <- .TreeForTaxa(tree, names(dataset))

  ret <- FastCharacterLength(tree, dataset)
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(dataset, "index")]
  }
}

#' Deprecated function
#' @keywords internal
#' @export
FitchSteps <- function(tree, dataset) {
  .Deprecated("CharacterLength")
  CharacterLength(tree, dataset, compress = TRUE)
}

#' @describeIn CharacterLength Do not perform checks.  Use with care: may cause
#' erroneous results or software crash if variables are in the incorrect format.
FastCharacterLength <- function(tree, dataset) {
  at <- attributes(dataset)
  if (is.null(at$levels) || ncol(at$contrast) == 0L) {
    return(rep(0L, at$nr))
  }
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  ts_char_steps(tree[["edge"]], at$contrast, tip_data,
                .ScaleWeight(at$weight), at$levels)
}

#' Score a tree from prepared data
#'
#' `TreeScore()` and `EdgeListScore()` compute the parsimony score of a tree
#' with the native C++ Fitch engine, using a dataset prepared by
#' [`PrepareData()`].
#' `EdgeListScore()` is the low-level scorer used as the default `TreeScorer`
#' by the custom-search functions; `TreeScore()` is a convenience wrapper that
#' takes a `phylo` tree.  Most users want the higher-level [`TreeLength()`].
#'
#' @template labelledTreeParam
#' @param dataset A `ParsimonyData` object from [`PrepareData()`].
#'
#' @return `TreeScore()` returns the parsimony score of `tree` under the
#' optimality criterion baked into `dataset` (equal weights, implied weights or
#' profile parsimony).
#'
#' @seealso [`PrepareData()`]; user-facing scoring with [`TreeLength()`].
#'
#' @family tree scoring
#' @author Martin R. Smith
#' @keywords internal
#' @importFrom TreeTools RenumberEdges RenumberTips
#' @export
TreeScore <- function(tree, dataset) {
  if (!is.ParsimonyData(dataset)) {
    stop("`dataset` must be a ParsimonyData object; see ?PrepareData.")
  }
  nTaxa <- dataset[["nTip"]]
  if (nTaxa != length(tree[["tip.label"]])) {
    stop("Number of taxa in dataset (", nTaxa,
         ") not equal to number of tips in tree")
  }
  tree <- RenumberTips(tree, dataset[["tip.label"]])
  el <- RenumberEdges(tree[["edge"]][, 1], tree[["edge"]][, 2])
  # Return:
  EdgeListScore(el[[1]], el[[2]], dataset)
}

#' @describeIn TreeScore Low-level scorer taking parent and child vectors; the
#'   default `TreeScorer` for [`TreeSearch()`], [`Ratchet()`] and
#'   [`Jackknife()`].  Scores via the native `ts_fitch_score()` kernel, which
#'   handles the Brazeau-Guillerme-Smith inapplicable algorithm correctly
#'   (including the `{1-}` case that MorphyLib mis-scored).
#' @inheritParams RearrangeEdges
#' @param inPostorder Logical: are the edges already in postorder?  If `FALSE`
#'   (the default) they are reordered before scoring.
#' @param \dots Unused; for interface compatibility with custom `TreeScorer`s.
#' @return `EdgeListScore()` returns the parsimony score (a numeric).
#' @importFrom TreeTools Preorder PostorderOrder
#' @export
EdgeListScore <- function(parent, child, dataset, inPostorder = FALSE, ...) {
  if (!inPostorder) {
    edgeList <- Preorder(cbind(parent, child))
    edgeList <- edgeList[PostorderOrder(edgeList), , drop = FALSE]
    parent <- edgeList[, 1]
    child <- edgeList[, 2]
  }
  # Return:
  ts_fitch_score(cbind(parent, child),
                 dataset[["contrast"]], dataset[["tip_data"]],
                 dataset[["weight"]], dataset[["levels"]],
                 min_steps = dataset[["min_steps"]],
                 concavity = dataset[["concavity"]],
                 infoAmounts = dataset[["info_amounts"]])
}

