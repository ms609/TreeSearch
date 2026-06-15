#' Simulate a dataset under parsimony
#'
#' Generates a morphological dataset under a strict parsimony model.
#' Characters are initialized at their minimum step count, then extra steps
#' are allocated one at a time. Each added step must increase the Fitch
#' parsimony score of the character by exactly one -- no "masked" or
#' "overprinted" steps are permitted.
#'
#' Back-mutations (e.g. \eqn{0 \to 1 \to 0}{0 -> 1 -> 0}) are allowed
#' when they genuinely add to the parsimony score.
#'
#' When `concavity` is finite (implied weighting), characters that already
#' carry more homoplasy are more likely to receive additional extra steps,
#' mirroring the mathematical relationship described by the
#' \eqn{k / (k + e)}{k/(k+e)} fit function.
#'
#' @param tree A \code{\link[ape:read.tree]{phylo}} object. If non-binary,
#'   resolved to binary with a warning. If unrooted, rooted internally at an
#'   arbitrary node. If no edge lengths are present, uniform lengths are used.
#' @param nChar Integer vector: `nChar[1]` gives the number of 2-state
#'   characters, `nChar[2]` the number of 3-state characters, and so on.
#' @param nExtraSteps Single integer: total extra steps distributed one at a
#'   time across all characters.
#' @param concavity Implied weighting concavity constant. `Inf` (default)
#'   gives equal weights (uniform character selection). A finite positive
#'   number _k_ gives implied weighting, with selection probability
#'   proportional to `(k + e) / k`. `"profile"` uses profile parsimony
#'   weighting: selection probability is proportional to the inverse of the
#'   phylogenetic information at the current step count, computed via
#'   [StepInformation()] after initialization.
#' @param rootState Integer vector: initial state at the root node for each
#'   character (default `0L`). If length 1, the same root state is used for
#'   all characters. If length `sum(nChar)`, each character gets its own root
#'   state. Each root state must be in `0:(k-1)` where _k_ is the number of
#'   states for that character.
#' @param missing Controls which cells are replaced with the ambiguous token
#'   `?`. Missing data is applied _after_ the complete simulation, so
#'   attributes such as `extra_steps` and `saturated` reflect the underlying
#'   complete dataset. Accepted forms:
#'   \describe{
#'     \item{**Scalar** (0--1)}{Flat rate: each cell is independently missing
#'       with this probability.}
#'     \item{**List** with `taxon` and/or `character` components}{Per-taxon
#'       and/or per-character rates. Each component is a numeric vector of
#'       probabilities (0--1). `taxon` should be named (matching tip labels)
#'       or length `n_tip`; `character` should be length `sum(nChar)`. Per-cell
#'       probability is `1 - (1 - p_taxon) * (1 - p_char)`.}
#'     \item{**Matrix** (n_tip x total_chars)}{Per-cell probabilities (0--1).
#'       Rows are taxa (named to match tip labels, or in tip order);
#'       columns are characters.}
#'   }
#'   Default `0` (no missing data).
#'
#' @return A `phyDat` object with characters ordered by number of states
#'   (2-state first, then 3-state, and so on). Additional attributes:
#'   \describe{
#'     \item{`saturated`}{Logical vector: can each character accept another
#'       step? Computed at return for all characters.}
#'     \item{`steps_exhausted`}{Logical vector: was each character discovered
#'       to be saturated during the step-placement loop (i.e., selected for
#'       an extra step but no legal edge found)?}
#'     \item{`extra_steps`}{Integer vector: number of extra steps placed on
#'       each character.}
#'   }
#'
#' @examples
#' tree <- TreeTools::BalancedTree(8)
#' dataset <- ParsSim(tree, nChar = c(20L), nExtraSteps = 10L)
#' TreeLength(tree, dataset)
#'
#' # Implied weighting: steps concentrate on fewer characters
#' dataset_iw <- ParsSim(tree, nChar = c(40L), nExtraSteps = 30L,
#'                       concavity = 3)
#' attr(dataset_iw, "extra_steps")
#'
#' # Profile parsimony weighting
#' dataset_pp <- ParsSim(tree, nChar = c(20L), nExtraSteps = 15L,
#'                       concavity = "profile")
#' attr(dataset_pp, "extra_steps")
#'
#' # 20% missing data injected post-hoc
#' dataset_missing <- ParsSim(tree, nChar = c(20L), nExtraSteps = 10L,
#'                            missing = 0.2)
#'
#' # Per-taxon missing rates (fragmentary taxa)
#' dataset_taxon <- ParsSim(tree, nChar = c(20L), nExtraSteps = 10L,
#'                          missing = list(taxon = c(t1 = 0.8, t2 = 0.5)))
#'
#' @references \insertCite{Goloboff2018}{TreeSearch}
#' \insertAllCited{}
#' @importFrom TreeTools MakeTreeBinary MatrixToPhyDat Postorder RootNode
#' @importFrom TreeTools RootTree
#' @family tree scoring
#' @export
ParsSim <- function(tree,
                    nChar = c(100L),
                    nExtraSteps = 0L,
                    concavity = Inf,
                    rootState = 0L,
                    missing = 0) {
  # --- Validate inputs -------------------------------------------------------
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylo object")
  }
  nChar <- as.integer(nChar)
  nExtraSteps <- as.integer(nExtraSteps)
  if (any(nChar < 0L)) stop("`nChar` values must be non-negative")
  total_chars <- sum(nChar)
  if (total_chars == 0L) stop("`nChar` must specify at least one character")
  if (length(nExtraSteps) != 1L || nExtraSteps < 0L) {
    stop("`nExtraSteps` must be a single non-negative integer")
  }
  missing_spec <- .pars_sim_validate_missing(missing)

  use_profile <- identical(concavity, "profile")
  if (!use_profile) {
    # `Inf` is valid (equal weights); reject NaN, -Inf and finite non-positive
    # values, which otherwise slip past the `is.finite()` test below and are
    # silently treated as equal weights.
    if (!is.numeric(concavity) || length(concavity) != 1L || is.nan(concavity) ||
        identical(concavity, -Inf) || (is.finite(concavity) && concavity <= 0)) {
      stop("`concavity` must be a positive number (or Inf for equal weights, ",
           "or \"profile\" for profile parsimony)")
    }
  }
  use_iw <- !use_profile && is.finite(concavity)

  # --- Prepare tree ----------------------------------------------------------
  tree_info <- .pars_sim_prepare_tree(tree)

  # --- Determine state counts per character ----------------------------------
  n_states_vec <- rep(seq_along(nChar) + 1L, times = nChar)

  # --- Validate and expand rootState ------------------------------------------
  rootState <- as.integer(rootState)
  if (length(rootState) == 1L) {
    rootState <- rep(rootState, total_chars)
  } else if (length(rootState) != total_chars) {
    stop("`rootState` must have length 1 or sum(nChar) (= ", total_chars, ")")
  }
  bad <- which(rootState < 0L | rootState >= n_states_vec)
  if (length(bad) > 0L) {
    stop("`rootState[", bad[1], "]` = ", rootState[bad[1]],
         " is out of range for a ", n_states_vec[bad[1]], "-state character",
         " (must be 0 to ", n_states_vec[bad[1]] - 1L, ")")
  }

  # --- Initialize characters -------------------------------------------------
  char_states <- vector("list", total_chars)
  char_scores <- integer(total_chars)
  for (i in seq_len(total_chars)) {
    init <- .pars_sim_init_char(tree_info, n_states_vec[i], rootState[i])
    char_states[[i]] <- init$node_states
    char_scores[i] <- init$score
  }

  # --- Compute info profiles for profile mode --------------------------------
  info_profiles <- NULL
  if (use_profile) {
    n_tip <- tree_info$n_tip
    info_profiles <- vector("list", total_chars)
    for (i in seq_len(total_chars)) {
      tip_states_i <- char_states[[i]][seq_len(n_tip)]
      info_profiles[[i]] <- StepInformation(tip_states_i)
    }
  }

  # --- Extra step loop -------------------------------------------------------
  extra_steps <- integer(total_chars)
  steps_exhausted <- logical(total_chars)

  if (nExtraSteps > 0L) {
    steps_placed <- 0L
    while (steps_placed < nExtraSteps) {
      available <- which(!steps_exhausted)
      if (length(available) == 0L) {
        warning("All characters saturated after ", steps_placed, " of ",
                nExtraSteps, " extra steps.")
        break
      }

      # Select character
      char_idx <- .pars_sim_select_char(available, extra_steps, concavity,
                                        use_iw, use_profile, char_scores,
                                        info_profiles)

      # Find legal edges
      legal <- .pars_sim_legal_edges(char_states[[char_idx]], tree_info,
                                     char_scores[char_idx],
                                     n_states_vec[char_idx])

      if (is.null(legal)) {
        steps_exhausted[char_idx] <- TRUE
        next
      }

      # Sample legal move weighted by edge length
      move_idx <- .safe_sample_idx(nrow(legal), prob = legal$edge_length)

      # Apply transition
      char_states[[char_idx]] <- .pars_sim_apply_transition(
        char_states[[char_idx]], tree_info,
        legal$edge_idx[move_idx], legal$target_state[move_idx]
      )
      char_scores[char_idx] <- char_scores[char_idx] + 1L
      extra_steps[char_idx] <- extra_steps[char_idx] + 1L
      steps_placed <- steps_placed + 1L

      # In profile mode, mark exhausted when info drops to 0
      if (use_profile) {
        profile <- info_profiles[[char_idx]]
        step_name <- as.character(char_scores[char_idx])
        if (!(step_name %in% names(profile)) ||
            profile[step_name] <= 0) {
          steps_exhausted[char_idx] <- TRUE
        }
      }
    }
  }

  # --- Build phyDat ----------------------------------------------------------
  n_tip <- tree_info$n_tip
  tip_matrix <- vapply(char_states, function(ns) ns[seq_len(n_tip)],
                       integer(n_tip))
  rownames(tip_matrix) <- tree_info$tip_labels

  prob_matrix <- .pars_sim_build_missing_matrix(
    missing_spec, n_tip, total_chars, tree_info$tip_labels
  )

  if (!is.null(prob_matrix)) {
    char_matrix <- matrix(as.character(tip_matrix), nrow = n_tip,
                          dimnames = dimnames(tip_matrix))
    is_missing <- matrix(runif(n_tip * total_chars), nrow = n_tip) < prob_matrix
    char_matrix[is_missing] <- "?"
    result <- MatrixToPhyDat(char_matrix)
  } else {
    result <- MatrixToPhyDat(tip_matrix)
  }

  # --- Calculate saturation for all characters --------------------------------
  saturated <- vapply(seq_len(total_chars), function(i) {
    is.null(.pars_sim_legal_edges(char_states[[i]], tree_info,
                                  char_scores[i], n_states_vec[i]))
  }, logical(1))

  attr(result, "saturated") <- saturated
  attr(result, "steps_exhausted") <- steps_exhausted
  attr(result, "extra_steps") <- extra_steps

  result
}


# --- Internal helpers --------------------------------------------------------

#' Prepare a tree for simulation
#' @return Named list: edge (postorder matrix), edge_length, n_tip, n_node,
#'   root, children (list of child-node vectors), tip_labels.
#' @keywords internal
#' @noRd
.pars_sim_prepare_tree <- function(tree) {
  if (!ape::is.rooted(tree)) {
    tree <- RootTree(tree, tree[["tip.label"]][1])
  }
  if (!ape::is.binary(tree)) {
    warning("Resolving non-binary tree to binary.")
    tree <- MakeTreeBinary(tree)
  }

  tree <- Postorder(tree)
  edge <- tree[["edge"]]
  n_tip <- length(tree[["tip.label"]])
  n_node <- n_tip + tree[["Nnode"]]

  edge_length <- tree[["edge.length"]]
  if (is.null(edge_length)) {
    edge_length <- rep(1, nrow(edge))
  }

  children <- vector("list", n_node)
  for (i in seq_len(n_node)) children[[i]] <- integer(0)
  for (i in seq_len(nrow(edge))) {
    p <- edge[i, 1]
    children[[p]] <- c(children[[p]], edge[i, 2])
  }

  list(
    edge = edge,
    edge_length = edge_length,
    n_tip = n_tip,
    n_node = n_node,
    root = RootNode(tree),
    children = children,
    tip_labels = tree[["tip.label"]]
  )
}


#' Fitch parsimony score for a single character
#'
#' Pure R Fitch downpass using bit-vector state sets.
#' @param tip_states Integer vector of states (0-indexed) for tips 1..n_tip.
#' @param tree_info List from `.pars_sim_prepare_tree()`.
#' @return Integer parsimony score.
#' @keywords internal
#' @noRd
.pars_sim_fitch_score <- function(tip_states, tree_info) {
  n_tip <- tree_info$n_tip
  n_node <- tree_info$n_node
  edge <- tree_info$edge

  sets <- integer(n_node)
  sets[seq_len(n_tip)] <- bitwShiftL(1L, tip_states[seq_len(n_tip)])

  score <- 0L
  for (i in seq_len(nrow(edge))) {
    p <- edge[i, 1]
    ch <- edge[i, 2]
    if (sets[p] == 0L) {
      sets[p] <- sets[ch]
    } else {
      inter <- bitwAnd(sets[p], sets[ch])
      if (inter > 0L) {
        sets[p] <- inter
      } else {
        sets[p] <- bitwOr(sets[p], sets[ch])
        score <- score + 1L
      }
    }
  }
  score
}


#' Initialize a character with minimum steps
#'
#' Sets all nodes to `root_state`, then places `n_states - 1` transitions on
#' randomly selected edges (weighted by edge length) to introduce each state.
#' @return List: `node_states` (integer vector, length n_node), `score`.
#' @keywords internal
#' @noRd
.pars_sim_init_char <- function(tree_info, n_states, root_state) {
  node_states <- rep(as.integer(root_state), tree_info$n_node)
  edge <- tree_info$edge

  other_states <- setdiff(seq.int(0L, n_states - 1L), root_state)
  for (new_state in other_states) {
    # Edges where both endpoints share the same state
    unmarked <- which(node_states[edge[, 1]] == node_states[edge[, 2]])
    weights <- tree_info$edge_length[unmarked]
    idx <- unmarked[.safe_sample_idx(length(unmarked), prob = weights)]

    node_states <- .pars_sim_apply_transition(node_states, tree_info, idx,
                                              new_state)
  }

  list(node_states = node_states, score = n_states - 1L)
}


#' Find contiguous region of same-state nodes below a start node
#'
#' DFS from `start_node` following edges where parent and child share the
#' same state.
#' @return List: `region` (all node indices), `tips` (tip-only indices),
#'   `boundary_states` (states of nodes just outside the region).
#' @keywords internal
#' @noRd
.pars_sim_find_region <- function(node_states, tree_info, start_node) {
  children <- tree_info$children
  n_tip <- tree_info$n_tip
  state <- node_states[start_node]

  region <- integer(0)
  tips <- integer(0)
  boundary_states <- integer(0)
  stack <- start_node

  while (length(stack) > 0L) {
    node <- stack[length(stack)]
    stack <- stack[-length(stack)]
    region <- c(region, node)
    if (node <= n_tip) {
      tips <- c(tips, node)
    } else {
      for (ch in children[[node]]) {
        if (node_states[ch] == state) {
          stack <- c(stack, ch)
        } else {
          boundary_states <- c(boundary_states, node_states[ch])
        }
      }
    }
  }

  list(region = region, tips = tips, boundary_states = boundary_states)
}


#' Find all legal (edge, target-state) moves for one character
#'
#' For each unmarked edge (endpoints share state), tries each possible
#' target state. Uses a boundary prefilter followed by Fitch verification.
#' @return Data frame with columns `edge_idx`, `target_state`, `edge_length`,
#'   or NULL if no legal moves.
#' @keywords internal
#' @noRd
.pars_sim_legal_edges <- function(node_states, tree_info, current_score,
                                  n_states) {
  edge <- tree_info$edge
  n_edge <- nrow(edge)

  edge_idx_out <- integer(0)
  target_state_out <- integer(0)
  edge_length_out <- numeric(0)
  all_states <- seq.int(0L, n_states - 1L)

  for (i in seq_len(n_edge)) {
    p <- edge[i, 1]
    ch <- edge[i, 2]

    # Only consider unmarked edges
    if (node_states[p] != node_states[ch]) next

    current_state <- node_states[ch]
    info <- .pars_sim_find_region(node_states, tree_info, ch)
    targets <- setdiff(all_states, current_state)

    for (t in targets) {
      # Boundary prefilter: if a boundary child already has target state,
      # the transition would eliminate an existing step → skip
      if (t %in% info$boundary_states) next

      # Fitch verify
      new_tip_states <- node_states[seq_len(tree_info$n_tip)]
      new_tip_states[info$tips] <- t
      new_score <- .pars_sim_fitch_score(new_tip_states, tree_info)

      if (new_score == current_score + 1L) {
        edge_idx_out <- c(edge_idx_out, i)
        target_state_out <- c(target_state_out, t)
        edge_length_out <- c(edge_length_out, tree_info$edge_length[i])
      }
    }
  }

  if (length(edge_idx_out) == 0L) return(NULL)

  data.frame(edge_idx = edge_idx_out,
             target_state = target_state_out,
             edge_length = edge_length_out)
}


#' Apply a transition on an edge
#'
#' Changes the child node and its contiguous same-state region to
#' `new_state`.
#' @return Updated `node_states` vector.
#' @keywords internal
#' @noRd
.pars_sim_apply_transition <- function(node_states, tree_info, edge_idx,
                                       new_state) {
  child_node <- tree_info$edge[edge_idx, 2]
  info <- .pars_sim_find_region(node_states, tree_info, child_node)
  node_states[info$region] <- new_state
  node_states
}


#' Select a character for the next extra step
#' @keywords internal
#' @noRd
.pars_sim_select_char <- function(available, extra_steps, concavity, use_iw,
                                  use_profile = FALSE, char_scores = NULL,
                                  info_profiles = NULL) {
  if (length(available) == 1L) return(available)

  if (use_profile) {
    # Weight ∝ 1 / info_amount at current step count
    weights <- vapply(available, function(i) {
      profile <- info_profiles[[i]]
      step_name <- as.character(char_scores[i])
      if (step_name %in% names(profile)) {
        info <- profile[step_name]
        if (info > 0) return(1.0 / info)
      }
      0
    }, double(1))
    # If all weights are 0, all available characters are info-saturated
    if (all(weights == 0)) return(available[1L])
    available[sample.int(length(available), 1L, prob = weights)]
  } else if (use_iw) {
    weights <- (concavity + extra_steps[available]) / concavity
    available[sample.int(length(available), 1L, prob = weights)]
  } else {
    available[sample.int(length(available), 1L)]
  }
}


#' Sample a single index, safe for length-1 vectors
#' @keywords internal
#' @noRd
.safe_sample_idx <- function(n, prob = NULL) {
  if (n == 1L) return(1L)
  if (!is.null(prob)) {
    # Edge lengths drive the weights; a tree with all-zero (or absent /
    # undefined) branch lengths leaves every candidate edge with weight 0,
    # for which sample.int() errors "too few positive probabilities". Treat
    # such edges as equiprobable instead.
    prob[is.na(prob)] <- 0
    if (!any(prob > 0)) {
      prob <- NULL
    }
  }
  sample.int(n, 1L, prob = prob)
}


#' Validate and parse the `missing` argument
#'
#' Returns a list with `type` ("none", "scalar", "list", "matrix") and
#' the parsed value.
#' @keywords internal
#' @noRd
.pars_sim_validate_missing <- function(missing) {
  if (is.matrix(missing)) {
    if (!is.numeric(missing)) stop("`missing` matrix must be numeric")
    if (any(is.na(missing)) || any(missing < 0) || any(missing > 1)) {
      stop("`missing` matrix values must be between 0 and 1")
    }
    return(list(type = "matrix", value = missing))
  }

  if (is.list(missing)) {
    valid_names <- c("taxon", "character")
    bad <- setdiff(names(missing), valid_names)
    if (length(bad) > 0L) {
      stop("`missing` list may only contain 'taxon' and/or 'character' ",
           "components; found: ", paste(bad, collapse = ", "))
    }
    if (length(missing) == 0L ||
        !any(valid_names %in% names(missing))) {
      stop("`missing` list must contain at least one of 'taxon' or 'character'")
    }
    for (comp in valid_names) {
      if (comp %in% names(missing)) {
        v <- missing[[comp]]
        if (!is.numeric(v) || any(is.na(v)) || any(v < 0) || any(v > 1)) {
          stop("`missing$", comp, "` must be a numeric vector with ",
               "values between 0 and 1")
        }
      }
    }
    return(list(type = "list", value = missing))
  }

  # Scalar case
  missing <- as.double(missing)
  if (length(missing) != 1L || is.na(missing) || missing < 0 || missing > 1) {
    stop("`missing` must be a number between 0 and 1, a list, or a matrix")
  }
  if (missing == 0) return(list(type = "none"))
  list(type = "scalar", value = missing)
}


#' Build a per-cell probability matrix from a missing specification
#'
#' @return A n_tip × total_chars matrix of probabilities, or NULL if no
#'   missing data should be applied.
#' @keywords internal
#' @noRd
.pars_sim_build_missing_matrix <- function(spec, n_tip, total_chars,
                                            tip_labels) {
  if (spec$type == "none") return(NULL)

  if (spec$type == "scalar") {
    return(matrix(spec$value, nrow = n_tip, ncol = total_chars))
  }

  if (spec$type == "matrix") {
    mat <- spec$value
    if (!is.null(rownames(mat))) {
      # Reorder rows to match tip_labels
      if (!all(tip_labels %in% rownames(mat))) {
        stop("`missing` matrix row names must include all tip labels")
      }
      mat <- mat[tip_labels, , drop = FALSE]
    }
    if (nrow(mat) != n_tip || ncol(mat) != total_chars) {
      stop("`missing` matrix must have ", n_tip, " rows (taxa) and ",
           total_chars, " columns (characters)")
    }
    return(mat)
  }

  # List case: combine taxon and character rates
  miss <- spec$value
  p_taxon <- rep(0, n_tip)
  if ("taxon" %in% names(miss)) {
    tv <- miss$taxon
    if (!is.null(names(tv))) {
      if (!all(names(tv) %in% tip_labels)) {
        stop("Names in `missing$taxon` must be valid tip labels")
      }
      # Named: match to tip labels; unlisted taxa get 0
      p_taxon[match(names(tv), tip_labels)] <- tv
    } else {
      if (length(tv) != n_tip) {
        stop("`missing$taxon` must be named or have length ", n_tip)
      }
      p_taxon <- tv
    }
  }

  p_char <- rep(0, total_chars)
  if ("character" %in% names(miss)) {
    cv <- miss$character
    if (length(cv) != total_chars) {
      stop("`missing$character` must have length ", total_chars)
    }
    p_char <- cv
  }

  # p_cell = 1 - (1 - p_taxon) * (1 - p_char)
  prob_mat <- 1 - outer(1 - p_taxon, 1 - p_char)

  if (all(prob_mat == 0)) return(NULL)
  prob_mat
}
