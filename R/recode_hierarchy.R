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
RecodeHierarchy <- function(dataset, hierarchy) {
  ValidateHierarchy(hierarchy, dataset)

  idx <- attr(dataset, "index")
  allLevels <- attr(dataset, "allLevels")
  nChar <- length(idx)
  nTip <- length(dataset)

  # Original character matrix (taxon × char), as token strings
  origMat <- do.call(rbind, lapply(dataset, function(x) {
    allLevels[x[idx]]
  }))

  .RecodeBlock <- function(node) {
    ctrl <- node$controlling
    deps <- node$dependents

    if (length(node$children) > 0L) {
      stop("Nested hierarchies not yet supported in RecodeHierarchy(). ",
           "Block controlled by character ", ctrl, " has sub-hierarchies.")
    }

    # Informative levels for each secondary (exclude "-" and "?")
    secLevels <- lapply(deps, function(d) {
      sort(setdiff(unique(origMat[, d]), c("-", "?")))
    })
    secNStates <- vapply(secLevels, length, integer(1))

    nPresent <- prod(secNStates)
    nStates <- nPresent + 1L
    nSec <- length(deps)

    if (nStates > 32L) {
      warning(sprintf(
        paste0("Hierarchy block controlled by character %d produces %d states ",
               "(> 32). Large state spaces may be slow."),
        ctrl, nStates
      ))
    }

    # All present-state combinations (expand.grid: first dim varies fastest)
    if (nSec > 0L) {
      comboGrid <- as.matrix(expand.grid(
        lapply(secLevels, seq_along)
      ))
    } else {
      # No secondaries: 2 states (absent + one present)
      comboGrid <- matrix(integer(0), nrow = 1L, ncol = 0L)
    }

    # --- Cost matrix ---
    gainCost <- nSec + 1L
    cm <- matrix(0, nStates, nStates)
    for (i in seq_len(nStates)) {
      for (j in seq_len(nStates)) {
        if (i == j) next
        if (i == 1L) {
          cm[i, j] <- gainCost  # absent → present
        } else if (j == 1L) {
          cm[i, j] <- 1         # present → absent
        } else {
          # Hamming distance between present combinations
          cm[i, j] <- sum(comboGrid[i - 1L, ] != comboGrid[j - 1L, ])
        }
      }
    }

    # --- Tip states ---
    tipStates <- integer(nTip)
    for (t in seq_len(nTip)) {
      pri <- origMat[t, ctrl]

      if (pri == "?") {
        tipStates[t] <- -1L  # fully ambiguous
        next
      }
      if (pri == "0" || pri == "-") {
        tipStates[t] <- 0L   # absent
        next
      }
      # Primary present: encode secondary combination
      if (nSec == 0L) {
        tipStates[t] <- 1L   # only present state
        next
      }

      secVals <- origMat[t, deps]
      anyUnknown <- FALSE
      levelIndices <- integer(nSec)

      for (s in seq_len(nSec)) {
        if (secVals[s] %in% c("-", "?")) {
          anyUnknown <- TRUE
          break
        }
        mi <- match(secVals[s], secLevels[[s]])
        if (is.na(mi)) {
          anyUnknown <- TRUE
          break
        }
        levelIndices[s] <- mi
      }

      if (anyUnknown) {
        tipStates[t] <- -2L  # present, unknown combination
        next
      }

      # Mixed-radix encoding (first dim varies fastest, matching expand.grid)
      rowIdx <- 1L
      multiplier <- 1L
      for (s in seq_len(nSec)) {
        rowIdx <- rowIdx + (levelIndices[s] - 1L) * multiplier
        multiplier <- multiplier * secNStates[s]
      }
      tipStates[t] <- rowIdx  # 1-based present state = Sankoff state index
    }

    list(
      n_states = nStates,
      cost_matrix = cm,
      tip_states = tipStates,
      forced_root_state = -1L,
      block_chars = c(ctrl, deps)
    )
  }

  blocks <- lapply(hierarchy, .RecodeBlock)

  hChars <- HierarchyChars(hierarchy)
  nonH <- setdiff(seq_len(nChar), hChars)

  list(
    sankoff_chars = blocks,
    non_hierarchy_indices = nonH
  )
}
