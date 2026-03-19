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
#' \insertCite{Goloboff2021b}{TreeSearch} approaches to inapplicable
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
#' @seealso [MaximizeParsimony()], [hierarchy_from_names()]
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
  controlling_indices <- suppressWarnings(as.integer(names(args)))
  if (anyNA(controlling_indices)) {
    stop("Controlling character names must be integer indices.")
  }

  lapply(seq_along(args), function(i) {
    ctrl <- controlling_indices[i]
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
        sub_ctrl <- suppressWarnings(as.integer(nms[j]))
        if (is.na(sub_ctrl)) {
          stop("Sub-hierarchy names must be integer character indices, got '",
               nms[j], "'.")
        }
        # The sub-controlling character is also a dependent of this block
        dependents <- c(dependents, sub_ctrl)
        children <- c(children, list(.ParseOneBlock(sub_ctrl, val[[j]])))
      }
    }
    return(list(
      controlling = as.integer(ctrl),
      dependents = dependents,
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
    leaf_deps <- setdiff(
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
#' @export
validate_hierarchy <- function(hierarchy, dataset) {
  if (!inherits(hierarchy, "CharacterHierarchy")) {
    stop("`hierarchy` must be a CharacterHierarchy object.")
  }
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }

  n_char <- length(attr(dataset, "index"))
  all_levels <- attr(dataset, "allLevels")
  levels <- attr(dataset, "levels")
  contrast <- attr(dataset, "contrast")

  # Identify the inapplicable token
  inapp_token <- "-"
  if (!inapp_token %in% all_levels) {
    stop("Dataset does not contain an inapplicable token ('-').")
  }

  # Build the original character matrix
  idx <- attr(dataset, "index")
  orig_mat <- do.call(rbind, lapply(dataset, function(x) {
    all_levels[x[idx]]
  }))

  # Identify the "0" state (absence) in the controlling primary
  absence_state <- "0"

  # Track all characters claimed by any block

  claimed <- integer(0)

  .ValidateBlock <- function(node, depth = 1L) {
    ctrl <- node$controlling
    deps <- node$dependents

    # Check indices exist
    all_idx <- c(ctrl, deps)
    bad <- all_idx[all_idx < 1L | all_idx > n_char]
    if (length(bad) > 0L) {
      stop(sprintf(
        "Character index(es) %s out of range [1, %d].",
        paste(bad, collapse = ", "), n_char
      ))
    }

    # Check no double-claiming
    overlap <- intersect(all_idx, claimed)
    if (length(overlap) > 0L) {
      stop(sprintf(
        "Character(s) %s appear in multiple hierarchy blocks.",
        paste(overlap, collapse = ", ")
      ))
    }
    claimed <<- c(claimed, all_idx)

    # Check controlling character is binary (has exactly states "0" and "1",
    # possibly with inapplicable/missing)
    ctrl_vals <- unique(orig_mat[, ctrl])
    ctrl_informative <- setdiff(ctrl_vals, c("?", "-"))
    if (!all(ctrl_informative %in% c("0", "1"))) {
      stop(sprintf(
        paste0("Controlling character %d must be binary (states '0' and '1'),",
               " but has states: %s."),
        ctrl, paste(ctrl_informative, collapse = ", ")
      ))
    }

    # Check secondaries are "-" where controlling is "0"
    absent_taxa <- which(orig_mat[, ctrl] == absence_state)
    if (length(absent_taxa) > 0L) {
      for (d in deps) {
        dep_vals <- orig_mat[absent_taxa, d]
        bad_taxa <- which(!dep_vals %in% c("-", "?"))
        if (length(bad_taxa) > 0L) {
          bad_names <- rownames(orig_mat)[absent_taxa[bad_taxa]]
          stop(sprintf(
            paste0("Secondary character %d has non-inapplicable values for ",
                   "taxa where controlling character %d is absent: %s."),
            d, ctrl, paste(head(bad_names, 5), collapse = ", ")
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
#' @param char_names Character vector of names, one per original character.
#'
#' @return A [`CharacterHierarchy`] object, or `NULL` if no hierarchy is
#'   detected.
#'
#' @examples
#' names <- c("sup_tail", "sub_tail_colour", "sub_tail_shape",
#'             "sup_wing", "sub_wing_venation", "eyes")
#' hierarchy_from_names(names)
#'
#' @family tree scoring
#' @seealso [CharacterHierarchy()]
#' @export
hierarchy_from_names <- function(char_names) {
  if (!is.character(char_names) || length(char_names) == 0L) {
    stop("`char_names` must be a non-empty character vector.")
  }

  # Find sup_ and sub_ characters
  sup_idx <- grep("^sup_", char_names)
  sub_idx <- grep("^sub_", char_names)

  if (length(sup_idx) == 0L) {
    return(NULL)
  }

  # Extract tags
  sup_tags <- sub("^sup_", "", char_names[sup_idx])
  sub_tags_full <- sub("^sub_", "", char_names[sub_idx])
  # The tag is the first component before any additional underscore-suffix
  # e.g. "sub_tail_colour" → tag = "tail"
  sub_tags <- sub("_.*", "", sub_tags_full)

  # Build mapping: tag → controlling index, tag → dependent indices
  tag_to_sup <- setNames(sup_idx, sup_tags)

  # Group sub characters by tag
  tag_to_subs <- split(sub_idx, sub_tags)

  # Check for sub_ characters referencing nonexistent sup_ tags
  orphan_tags <- setdiff(names(tag_to_subs), sup_tags)
  if (length(orphan_tags) > 0L) {
    warning(sprintf(
      "sub_ characters reference tags with no corresponding sup_: %s",
      paste(orphan_tags, collapse = ", ")
    ))
  }

  # Detect nested hierarchies: a sub_ character that is also a sup_
  # Find sub_ chars that are also in sup_idx
  sub_also_sup <- intersect(sub_idx, sup_idx)

  # Build hierarchy
  # First pass: create flat blocks for all sup_ tags
  args <- list()
  for (tag in sup_tags) {
    ctrl <- tag_to_sup[[tag]]
    subs <- tag_to_subs[[tag]]
    if (is.null(subs)) subs <- integer(0)

    # Check which subs are themselves controlling (nested hierarchy)
    nested_subs <- intersect(subs, sup_idx)
    flat_subs <- setdiff(subs, sup_idx)

    if (length(nested_subs) == 0L) {
      # Simple block
      args[[as.character(ctrl)]] <- as.integer(subs)
    } else {
      # Nested: build list with named sub-hierarchies
      block <- as.list(as.integer(flat_subs))
      for (ns in nested_subs) {
        ns_tag <- sup_tags[sup_idx == ns]
        ns_subs <- tag_to_subs[[ns_tag]]
        if (is.null(ns_subs)) ns_subs <- integer(0)
        block[[as.character(ns)]] <- as.integer(ns_subs)
      }
      args[[as.character(ctrl)]] <- block
    }
  }

  # Filter out sup_ chars whose index also appears in sub_idx
  # (they'll be included as children of their parent)
  top_level_sup <- setdiff(sup_idx, sub_idx)
  if (length(top_level_sup) == 0L) {
    # All sup_ characters are also sub_ — circular or all nested.
    # Fall back to treating all as top-level with a warning.
    warning("All sup_ characters are also sub_ characters. ",
            "Treating all as top-level.")
    top_level_sup <- sup_idx
  }
  top_level_ctrls <- as.character(top_level_sup)
  args <- args[top_level_ctrls]

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
hierarchy_chars <- function(hierarchy) {
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
hierarchy_controlling <- function(hierarchy) {
  vapply(hierarchy, `[[`, integer(1), "controlling")
}


#' Build tip_labels matrix for HSJ scoring
#'
#' Converts a `phyDat` dataset into an integer matrix of per-tip per-character
#' state labels (0-based) for the HSJ C++ scoring function.
#'
#' @param dataset A `phyDat` object.
#' @return An integer matrix with `length(dataset)` rows (tips) and
#'   `length(attr(dataset, "index"))` columns (original characters).
#'   Each entry is a 0-based token index.
#' @keywords internal
#' @export
build_tip_labels <- function(dataset) {
  idx <- attr(dataset, "index")
  n_tip <- length(dataset)
  n_char <- length(idx)

  # dataset is a list of integer vectors (pattern indices per tip)
  # Expand via index to original characters, convert to 0-based
  mat <- matrix(0L, nrow = n_tip, ncol = n_char)
  for (t in seq_len(n_tip)) {
    pattern_tokens <- dataset[[t]]    # token indices for each pattern
    mat[t, ] <- pattern_tokens[idx] - 1L  # 0-based
  }
  mat
}


#' Convert CharacterHierarchy to list for C++
#'
#' Converts a [`CharacterHierarchy`] object into a flat list of hierarchy
#' blocks that can be passed to the C++ `ts_hsj_score()` bridge function.
#' Each block is a list with `primary` (0-based) and `secondaries` (0-based).
#'
#' @param hierarchy A [`CharacterHierarchy`] object.
#' @return A list of lists, each with elements `primary` (integer, 0-based)
#'   and `secondaries` (integer vector, 0-based).
#' @keywords internal
#' @export
hierarchy_to_blocks <- function(hierarchy) {
  .flatten_block <- function(node) {
    block <- list(
      primary = node$controlling - 1L,
      secondaries = node$dependents - 1L
    )
    child_blocks <- lapply(node$children, .flatten_block)
    c(list(block), unlist(child_blocks, recursive = FALSE))
  }
  unlist(lapply(hierarchy, .flatten_block), recursive = FALSE)
}


#' Compute non-hierarchy pattern weights
#'
#' Given a `phyDat` dataset and a [`CharacterHierarchy`], returns a weight
#' vector with hierarchy characters' contributions subtracted.
#' Patterns that appear only in hierarchy characters will have weight 0.
#'
#' @param dataset A `phyDat` object.
#' @param hierarchy A [`CharacterHierarchy`] object.
#'
#' @return An integer vector of adjusted pattern weights (same length as
#'   `attr(dataset, "weight")`).
#'
#' @keywords internal
#' @export
non_hierarchy_weights <- function(dataset, hierarchy) {
  w <- attr(dataset, "weight")
  idx <- attr(dataset, "index")
  h_chars <- hierarchy_chars(hierarchy)

  adjusted <- as.integer(w)
  for (ci in h_chars) {
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
#   non_hierarchy_weights: pattern weights for Fitch scoring (non-hierarchy
#     chars only, reflecting which free chars were sampled)
#   block_counts: integer vector (length = number of top-level blocks)
#     giving how many times each block was sampled (0/1 for jackknife,
#     0+ for bootstrap)
.HierarchicalResampleWeights <- function(dataset, hierarchy, bootstrap,
                                         proportion) {
  idx <- attr(dataset, "index")
  n_patterns <- length(attr(dataset, "weight"))
  n_chars <- length(idx)

  # Collect chars per top-level block (includes nested dependents)
  .CollectAll <- function(node) {
    c(node$controlling, node$dependents,
      unlist(lapply(node$children, .CollectAll)))
  }
  n_blocks <- length(hierarchy)
  block_chars <- lapply(hierarchy, function(node) unique(.CollectAll(node)))
  h_chars_set <- unique(unlist(block_chars))

  free_chars <- setdiff(seq_len(n_chars), h_chars_set)
  n_free <- length(free_chars)
  n_units <- n_free + n_blocks

  if (n_units < 2L) {
    # Degenerate: can't jackknife with < 2 units
    return(list(
      non_hierarchy_weights = non_hierarchy_weights(dataset, hierarchy),
      block_counts = rep(1L, n_blocks)
    ))
  }

  if (bootstrap) {
    sampled <- sample.int(n_units, n_units, replace = TRUE)
  } else {
    n_keep <- max(1L, ceiling(proportion * n_units))
    n_keep <- min(n_keep, n_units - 1L)
    sampled <- sample.int(n_units, n_keep, replace = FALSE)
  }

  unit_counts <- tabulate(sampled, nbins = n_units)

  # Non-hierarchy pattern weights from retained free chars
  nh_weights <- integer(n_patterns)
  for (i in seq_len(n_free)) {
    if (unit_counts[i] > 0L) {
      pat <- idx[free_chars[i]]
      nh_weights[pat] <- nh_weights[pat] + unit_counts[i]
    }
  }

  block_counts <- unit_counts[n_free + seq_len(n_blocks)]

  list(
    non_hierarchy_weights = nh_weights,
    block_counts = block_counts
  )
}
