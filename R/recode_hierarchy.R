#' Recode hierarchical characters as step-matrix characters
#'
#' Implements the x-transformation recoding of
#' \insertCite{Goloboff2021;textual}{TreeSearch}.
#' Each hierarchy block (one controlling primary character plus \eqn{n}
#' secondary characters) is combined into a single step-matrix character
#' with \eqn{\prod k_i + 1} states and an asymmetric cost matrix.
#'
#' @details
#' ## State encoding
#'
#' State 0 represents "primary absent".
#' States \eqn{1 \ldots \prod k_i} represent all possible combinations of
#' secondary character states (where \eqn{k_i} is the number of informative
#' states of secondary character \eqn{i}).
#'
#' ## Cost matrix
#'
#' - **Absent → present (gain):** cost = \eqn{n + 1}, where \eqn{n} is the
#'   number of secondary characters.
#' - **Present → absent (loss):** cost = 1.
#' - **Present → present:** Hamming distance (number of secondaries with
#'   different states).
#'
#' @param dataset A [`phyDat`][phangorn::phyDat] object.
#' @param hierarchy A [`CharacterHierarchy`] object.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`sankoff_chars`}{A list of per-block lists, each containing:
#'     \describe{
#'       \item{`n_states`}{Integer, number of states (absent + present combos).}
#'       \item{`cost_matrix`}{Numeric matrix (\code{n_states × n_states}),
#'         row-major: \code{cost_matrix[from, to]}.}
#'       \item{`tip_states`}{Integer vector (length \code{n_tip}, 0-based).
#'         0 = absent, 1..n_present = present combination,
#'         -1 = fully ambiguous (all states possible),
#'         -2 = present but unknown combination.}
#'       \item{`forced_root_state`}{Integer: -1 (unconstrained).}
#'       \item{`block_chars`}{Integer vector of original character indices
#'         (1-based) belonging to this block.}
#'     }
#'   }
#'   \item{`non_hierarchy_indices`}{Integer vector of original character
#'     indices (1-based) not in any hierarchy block.}
#' }
#'
#' @references
#' \insertAllCited{}
#' @family tree scoring
#' @seealso [CharacterHierarchy()], [MaximizeParsimony()]
#' @keywords internal
#' @export
recode_hierarchy <- function(dataset, hierarchy) {
  validate_hierarchy(hierarchy, dataset)

  idx <- attr(dataset, "index")
  all_levels <- attr(dataset, "allLevels")
  n_char <- length(idx)
  n_tip <- length(dataset)

  # Original character matrix (taxon × char), as token strings
  orig_mat <- do.call(rbind, lapply(dataset, function(x) {
    all_levels[x[idx]]
  }))

  .recode_block <- function(node) {
    ctrl <- node$controlling
    deps <- node$dependents

    if (length(node$children) > 0L) {
      stop("Nested hierarchies not yet supported in recode_hierarchy(). ",
           "Block controlled by character ", ctrl, " has sub-hierarchies.")
    }

    # Informative levels for each secondary (exclude "-" and "?")
    sec_levels <- lapply(deps, function(d) {
      sort(setdiff(unique(orig_mat[, d]), c("-", "?")))
    })
    sec_nstates <- vapply(sec_levels, length, integer(1))

    n_present <- prod(sec_nstates)
    n_states <- n_present + 1L
    n_sec <- length(deps)

    if (n_states > 32L) {
      warning(sprintf(
        paste0("Hierarchy block controlled by character %d produces %d states ",
               "(> 32). Large state spaces may be slow."),
        ctrl, n_states
      ))
    }

    # All present-state combinations (expand.grid: first dim varies fastest)
    if (n_sec > 0L) {
      combo_grid <- as.matrix(expand.grid(
        lapply(sec_levels, seq_along)
      ))
    } else {
      # No secondaries: 2 states (absent + one present)
      combo_grid <- matrix(integer(0), nrow = 1L, ncol = 0L)
    }

    # --- Cost matrix ---
    gain_cost <- n_sec + 1L
    cm <- matrix(0, n_states, n_states)
    for (i in seq_len(n_states)) {
      for (j in seq_len(n_states)) {
        if (i == j) next
        if (i == 1L) {
          cm[i, j] <- gain_cost  # absent → present
        } else if (j == 1L) {
          cm[i, j] <- 1          # present → absent
        } else {
          # Hamming distance between present combinations
          cm[i, j] <- sum(combo_grid[i - 1L, ] != combo_grid[j - 1L, ])
        }
      }
    }

    # --- Tip states ---
    tip_states <- integer(n_tip)
    for (t in seq_len(n_tip)) {
      pri <- orig_mat[t, ctrl]

      if (pri == "?") {
        tip_states[t] <- -1L  # fully ambiguous
        next
      }
      if (pri == "0" || pri == "-") {
        tip_states[t] <- 0L   # absent
        next
      }
      # Primary present: encode secondary combination
      if (n_sec == 0L) {
        tip_states[t] <- 1L   # only present state
        next
      }

      sec_vals <- orig_mat[t, deps]
      any_unknown <- FALSE
      level_indices <- integer(n_sec)

      for (s in seq_len(n_sec)) {
        if (sec_vals[s] %in% c("-", "?")) {
          any_unknown <- TRUE
          break
        }
        mi <- match(sec_vals[s], sec_levels[[s]])
        if (is.na(mi)) {
          any_unknown <- TRUE
          break
        }
        level_indices[s] <- mi
      }

      if (any_unknown) {
        tip_states[t] <- -2L  # present, unknown combination
        next
      }

      # Mixed-radix encoding (first dim varies fastest, matching expand.grid)
      row_idx <- 1L
      multiplier <- 1L
      for (s in seq_len(n_sec)) {
        row_idx <- row_idx + (level_indices[s] - 1L) * multiplier
        multiplier <- multiplier * sec_nstates[s]
      }
      tip_states[t] <- row_idx  # 1-based present state = Sankoff state index
    }

    list(
      n_states = n_states,
      cost_matrix = cm,
      tip_states = tip_states,
      forced_root_state = -1L,
      block_chars = c(ctrl, deps)
    )
  }

  blocks <- lapply(hierarchy, .recode_block)

  h_chars <- hierarchy_chars(hierarchy)
  non_h <- setdiff(seq_len(n_char), h_chars)

  list(
    sankoff_chars = blocks,
    non_hierarchy_indices = non_h
  )
}
