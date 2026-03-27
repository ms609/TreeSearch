#' Score-spectrum coverage estimate for parsimony search
#'
#' `ScoreSpectrum()` applies Chao1-style abundance-based richness estimation
#' \insertCite{Chao1984}{TreeSearch} to the distribution of per-replicate
#' parsimony scores returned by [MaximizeParsimony()].  Treating each distinct
#' score value as a "species" and the number of replicates that found it as its
#' "abundance", the estimator quantifies how thoroughly the search has explored
#' the parsimony landscape.
#'
#' The **sample coverage** (Good-Turing estimator)
#' \insertCite{Good1953,Chao2012}{TreeSearch} is:
#' \deqn{\hat{C} = 1 - f_1 / n}
#' where \eqn{f_1} is the number of score levels seen exactly once and \eqn{n}
#' is the total number of replicates.  A coverage close to 1 indicates that
#' most of the accessible score landscape has been sampled; low coverage
#' suggests meaningful unexplored territory remains.
#'
#' The **Chao1 lower bound** on total score-level richness is:
#' \deqn{\hat{S} = S_{\mathrm{obs}} + \frac{f_1^2}{2 f_2}}
#' When \eqn{f_2 = 0} (no doubleton scores) the bias-corrected form
#' \eqn{f_1(f_1 - 1)/2} is used instead.
#'
#' @param trees A `multiPhylo` object returned by [MaximizeParsimony()], which
#'   must carry a `replicate_scores` attribute.  Alternatively, a numeric
#'   vector of per-replicate scores.
#' @param tol Numeric tolerance for binning floating-point scores.  Scores
#'   that differ by less than `tol` are treated as equal.  The default
#'   (`1e-4`) is suitable for implied-weights and profile-parsimony scores;
#'   use `0` for strict equality when working with equal-weights (integer)
#'   scores.
#'
#' @return A list of class `"ScoreSpectrum"` with components:
#'   \describe{
#'     \item{`n_replicates`}{Total completed replicates.}
#'     \item{`observed_levels`}{Distinct score values observed (\eqn{S_\mathrm{obs}}).}
#'     \item{`estimated_levels`}{Chao1 lower-bound estimate of total score
#'       levels (\eqn{\hat{S}}).}
#'     \item{`coverage`}{Good-Turing sample coverage (\eqn{\hat{C}}).}
#'     \item{`unseen_fraction`}{Estimated fraction of score levels not yet
#'       seen: \eqn{1 - S_\mathrm{obs}/\hat{S}}.}
#'     \item{`best_score`}{The lowest (best) score found.}
#'     \item{`best_score_reps`}{Number of replicates that reached the best
#'       score.}
#'     \item{`f`}{Named integer vector: \eqn{f_k} = number of score levels
#'       seen exactly \eqn{k} times (frequency spectrum).}
#'     \item{`replicate_scores`}{The raw per-replicate scores.}
#'   }
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library("TreeTools", quietly = TRUE)
#' data("Lobo", package = "TreeSearch")
#' result <- MaximizeParsimony(Lobo.phy, maxReplicates = 4L)
#' sp <- ScoreSpectrum(result)
#' print(sp)
#'
#' @family search utilities
#' @export
ScoreSpectrum <- function(trees, tol = 1e-4) {
  # Accept either a multiPhylo with attribute or a raw numeric vector
  if (inherits(trees, "multiPhylo")) {
    scores <- attr(trees, "replicate_scores")
    if (is.null(scores)) {
      stop("`trees` has no `replicate_scores` attribute. ",
           "Re-run MaximizeParsimony() with this version of TreeSearch.")
    }
  } else if (is.numeric(trees)) {
    scores <- trees
  } else {
    stop("`trees` must be a `multiPhylo` from MaximizeParsimony() or a ",
         "numeric vector of per-replicate scores.")
  }

  scores <- scores[is.finite(scores)]
  n <- length(scores)

  if (n < 2L) {
    return(structure(
      list(
        n_replicates = n,
        observed_levels = if (n == 0L) 0L else 1L,
        estimated_levels = NA_real_,
        coverage = NA_real_,
        unseen_fraction = NA_real_,
        best_score = if (n > 0L) min(scores) else NA_real_,
        best_score_reps = if (n > 0L) sum(scores == min(scores)) else 0L,
        f = integer(0L),
        replicate_scores = scores
      ),
      class = "ScoreSpectrum"
    ))
  }

  # Bin scores to handle floating-point equality
  if (tol > 0) {
    scores_binned <- round(scores / tol) * tol
  } else {
    scores_binned <- scores
  }

  # Frequency of each distinct score value (abundance vector)
  abundance <- tabulate(factor(scores_binned))
  s_obs <- length(abundance)   # distinct score levels observed

  # Frequency spectrum: f_k = number of score levels seen exactly k times
  max_k <- max(abundance)
  f_k <- tabulate(abundance, nbins = max_k)
  f1 <- if (max_k >= 1L) f_k[1L] else 0L
  f2 <- if (max_k >= 2L) f_k[2L] else 0L

  # Good-Turing sample coverage
  coverage <- 1.0 - f1 / n

  # Chao1 lower-bound estimate of total richness
  if (f2 > 0L) {
    s_hat <- s_obs + f1^2 / (2 * f2)
  } else if (f1 > 1L) {
    # Bias-corrected form when no doubletons
    s_hat <- s_obs + f1 * (f1 - 1L) / 2
  } else {
    # All observed levels are well-represented
    s_hat <- s_obs
  }

  unseen_fraction <- if (s_hat > 0) 1 - s_obs / s_hat else 0

  best_score <- min(scores_binned)
  best_score_reps <- sum(scores_binned <= best_score + tol)

  # Trim trailing zeros from f_k for a compact spectrum
  last_nonzero <- max(which(f_k > 0L), 0L)
  f_k_trimmed <- f_k[seq_len(last_nonzero)]
  names(f_k_trimmed) <- seq_len(last_nonzero)

  structure(
    list(
      n_replicates = n,
      observed_levels = s_obs,
      estimated_levels = s_hat,
      coverage = coverage,
      unseen_fraction = unseen_fraction,
      best_score = best_score,
      best_score_reps = best_score_reps,
      f = f_k_trimmed,
      replicate_scores = scores
    ),
    class = "ScoreSpectrum"
  )
}

#' @export
print.ScoreSpectrum <- function(x, ...) {
  if (is.na(x$coverage)) {
    cat("ScoreSpectrum: insufficient replicates (n =", x$n_replicates, ")\n")
    return(invisible(x))
  }
  cat(sprintf(
    "Score-spectrum coverage (n = %d replicates)\n",
    x$n_replicates
  ))
  cat(sprintf(
    "  Best score:         %.4g  (%d replicates)\n",
    x$best_score, x$best_score_reps
  ))
  cat(sprintf(
    "  Score levels seen:  %d  (est. total: %.1f)\n",
    x$observed_levels, x$estimated_levels
  ))
  cat(sprintf(
    "  Landscape coverage: %.1f%%",
    100 * x$coverage
  ))
  if (x$unseen_fraction > 0.01) {
    cat(sprintf("  (~%.0f%% of score levels unseen)", 100 * x$unseen_fraction))
  }
  cat("\n")
  if (length(x$f) > 0L) {
    cat("  Frequency spectrum (f_k): ")
    cat(paste0("f", names(x$f), "=", x$f, collapse = ", "), "\n")
  }
  invisible(x)
}
