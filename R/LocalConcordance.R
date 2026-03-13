#' Local concordance of an alignment
#'
#' `LocalConcordance()` computes an information-theoretic measure of regional
#' character agreement across alignment sites.  For each site, it reports the
#' Gaussian-weighted mean normalized mutual information (NMI) with neighbouring
#' sites, evaluated at multiple spatial scales.  Low scores flag regions where
#' characters are mutually inconsistent, analogous to the "garbage regions"
#' identified by PhyIN \insertCite{Maddison2024}{TreeSearch}.
#'
#' @param dataset A `phyDat` object.
#' @param sigma Numeric vector of spatial scales (in alignment columns) at
#'   which to compute local concordance.  Analogous to PhyIN's `-d` parameter
#'   (maximum neighbour distance), but continuous: each `sigma` value is the
#'   standard deviation of a Gaussian kernel.  The weight assigned to a
#'   neighbour at distance \eqn{d} is proportional to
#'   \eqn{\exp(-d^2 / 2\sigma^2)}, so the weight falls to
#'   \eqn{e^{-1/2} \approx 61\%} of its \eqn{d = 1} value at
#'   \eqn{d = \sigma}.  PhyIN's default of \code{-d 2} (uniform weight up to
#'   distance 2, zero beyond) is most closely approximated by
#'   \eqn{\sigma \approx 1.5}, at which the \eqn{d = 2} weight is roughly
#'   half the \eqn{d = 1} weight and negligible beyond \eqn{d = 3}.
#' @param truncate Positive number; the Gaussian window is truncated at
#'   `ceiling(truncate * max(sigma))` columns from the focal site.
#' @param normalize Logical; if `TRUE` (default), NMI scores are corrected for
#'   expected mutual information under statistical independence, using the
#'   analytical formula of \insertCite{Vinh2010;textual}{TreeSearch}.  Scores
#'   of 0 correspond to the expectation under independence; negative scores
#'   indicate less concordance than expected by chance.
#' @param internal_gaps Logical (default `TRUE`); analogous to PhyIN's `-e`
#'   flag.  If `TRUE`, gaps that are flanked on both sides by observed
#'   nucleotides within the same sequence (internal gaps, likely representing
#'   real insertion/deletion events) are recoded as a distinct extra state
#'   before computing NMI.  Terminal gaps (contiguous with the sequence end)
#'   are still treated as missing data.  When working on a concatenated
#'   alignment, the distinction between internal and terminal gaps is lost
#'   across locus boundaries; in that case consider `internal_gaps = FALSE`.
#' @param block_size Integer or `NULL` (default `NULL`).  Analogous to PhyIN's
#'   `-b` parameter.  If an integer, a second-pass regional scan is performed
#'   after computing per-site scores: a rolling window of `block_size` sites
#'   is slid across the alignment, and any site whose window contains a
#'   proportion of "noisy" sites \eqn{\ge} `threshold` is flagged for
#'   trimming.  The flagged sites are returned as the `"trim"` logical
#'   attribute of the result.  If `NULL`, the block scan is skipped.
#' @param threshold Proportion in \eqn{[0,\,1]} (default `0.5`); analogous to
#'   PhyIN's `-p` parameter.  The minimum proportion of "noisy" sites within a
#'   rolling block required to flag that block for trimming.  Increasing
#'   `threshold` makes the trimming more permissive.
#' @param noise_level Numeric (default `0`).  A site is classified as "noisy"
#'   for the block scan if its NMI score (at the smallest `sigma`) is
#'   \eqn{\le} `noise_level`.  The default of 0 is natural when
#'   `normalize = TRUE`, where 0 corresponds to the random-association
#'   baseline.
#'
#' @return A numeric matrix with one row per alignment column and one column
#'   per element of `sigma`.  Row names are column indices (as characters);
#'   column names are the `sigma` values.  The attribute `"pairNMI"` contains
#'   the pairwise NMI matrix over unique site patterns (dimensions
#'   `nr` \eqn{\times} `nr`).  If `block_size` is not `NULL`, the logical
#'   attribute `"trim"` indicates which sites the block scan recommends for
#'   deletion.
#'
#' @details
#' The score for site \eqn{i} at scale \eqn{\sigma} is:
#' \deqn{\mathrm{score}_\sigma(i) =
#'   \frac{\sum_{j \neq i} w_{ij} \, h^*_{ij} \, \mathrm{NMI}_{ij}}
#'        {\sum_{j \neq i} w_{ij} \, h^*_{ij}}}
#' where \eqn{w_{ij} = \phi(|j-i|;\, 0,\, \sigma)} is a Gaussian kernel
#' weight, \eqn{h^*_{ij} = \min\bigl(H(s_i), H(s_j)\bigr)} is the maximum
#' possible mutual information for the pair (used as an importance weight), and
#' \eqn{\mathrm{NMI}_{ij}} is the normalized mutual information between sites
#' \eqn{i} and \eqn{j}.  Sites whose neighbour window contains no informative
#' pairs receive `NA`.
#'
#' No tree is needed; concordance is measured entirely between characters.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' dataset <- congreveLamsdellMatrices[[1]]
#' lc <- LocalConcordance(dataset, sigma = c(1, 2, 5), internal_gaps = FALSE)
#' head(lc)
#' @seealso [PlotLocalConcordance()] for visualization.
#' @template MRS
#' @importFrom stats dnorm filter
#' @importFrom TreeDist entropy_int
#' @family split support functions
#' @export
LocalConcordance <- function(
  dataset,
  sigma         = c(0.5, 0.8, 1, 1.2, 1.5, 2, 3, 5),
  truncate      = 3,
  normalize     = TRUE,
  internal_gaps = TRUE,
  block_size    = NULL,
  threshold     = 0.5,
  noise_level   = 0
) {
  if (is.null(dataset)) {
    warning("Cannot compute local concordance without `dataset`.")
    return(NULL)
  }
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }

  # ---- 1. Extract raw attributes ----
  at         <- attributes(dataset)
  cont       <- at[["contrast"]]
  allLevels  <- at[["allLevels"]]
  idx_orig   <- at[["index"]]   # length L; values in 1..nr
  L          <- length(idx_orig)
  nr_orig    <- at[["nr"]]

  # Identify which tokens are gaps ("-") and which are ambiguous
  if ("-" %in% colnames(cont)) {
    cont[cont[, "-"] > 0, ] <- 1
  }
  ambiguous <- rowSums(cont) != 1   # logical, length = number of phyDat tokens
  gap_token <- if ("-" %in% allLevels) which(allLevels == "-") else integer(0)

  # Raw integer matrix: nTaxa × nr_orig (unique patterns)
  mat_raw <- matrix(as.integer(unlist(dataset)), length(dataset), byrow = TRUE)
  nTaxa   <- nrow(mat_raw)

  # ---- 2. Internal gap recoding (-e analogue) ----
  # Must work on the full L-column matrix to determine per-taxon flanking.
  if (internal_gaps && length(gap_token) > 0) {
    # Expand to full alignment: nTaxa × L
    full_mat <- mat_raw[, idx_orig, drop = FALSE]

    internalGapToken <- max(mat_raw, na.rm = TRUE) + 1L

    for (t in seq_len(nTaxa)) {
      row    <- full_mat[t, ]
      is_gap <- row %in% gap_token
      if (!any(is_gap) || all(is_gap)) next

      non_gap_pos <- which(!is_gap)
      first_ng    <- non_gap_pos[1L]
      last_ng     <- non_gap_pos[length(non_gap_pos)]

      # Internal: gap that lies strictly between first and last non-gap
      is_internal <- is_gap &
        (seq_len(L) > first_ng) &
        (seq_len(L) < last_ng)
      full_mat[t, is_internal] <- internalGapToken
    }

    # Re-derive unique patterns from the recoded full matrix
    colKeys    <- apply(full_mat, 2, paste, collapse = "\r")
    first_occ  <- which(!duplicated(colKeys))
    mat        <- full_mat[, first_occ, drop = FALSE]  # nTaxa × new_nr
    idx        <- match(colKeys, colKeys[first_occ])
    nr         <- length(first_occ)
    maxToken   <- internalGapToken
  } else {
    mat      <- mat_raw
    idx      <- idx_orig
    nr       <- nr_orig
    maxToken <- max(mat_raw, na.rm = TRUE)
    internalGapToken <- NA_integer_
  }

  # ---- 3. Standard preprocessing ----
  # Ambiguous tokens → NA  (internal gap token survives, as intended)
  orig_ambig <- which(ambiguous)   # only covers original tokens
  mat[mat %in% orig_ambig] <- NA_integer_

  # Singleton (autapomorphic) states → NA.
  # Exception: the internal gap token is NEVER dropped as a singleton — an
  # indel event in a single lineage is a real evolutionary signal.
  mat <- apply(mat, 2, function(x) {
    singletons <- which(tabulate(x, maxToken) == 1L)
    if (!is.na(internalGapToken)) {
      singletons <- setdiff(singletons, internalGapToken)
    }
    x[x %in% singletons] <- NA_integer_
    x
  })
  # mat: nTaxa × nr

  # ---- 4. Per-pattern entropy ----
  patH   <- double(nr)
  patTbl <- vector("list", nr)
  for (p in seq_len(nr)) {
    obs <- mat[!is.na(mat[, p]), p]
    if (length(obs) >= 2L) {
      tbl         <- tabulate(obs, maxToken)
      patH[p]     <- entropy_int(tbl)
      patTbl[[p]] <- tbl
    }
  }

  # ---- 5. Pairwise NMI matrix (nr × nr) ----
  nmi_mat   <- matrix(NA_real_, nr, nr)
  hBest_mat <- matrix(0,        nr, nr)

  for (pi in seq_len(nr)) {
    if (patH[pi] == 0) next
    for (pj in seq_len(pi - 1L)) {
      if (patH[pj] == 0) next

      aIJ <- !is.na(mat[, pi]) & !is.na(mat[, pj])
      cI  <- mat[aIJ, pi]
      cJ  <- mat[aIJ, pj]
      if (length(cI) < 2L) next

      tblI  <- tabulate(cI, maxToken)
      tblJ  <- tabulate(cJ, maxToken)
      hI    <- entropy_int(tblI)
      hJ    <- entropy_int(tblJ)
      hBest <- min(hI, hJ)
      if (hBest < sqrt(.Machine$double.eps)) next

      joint <- tabulate(cJ + (cI - 1L) * maxToken, maxToken * maxToken)
      hIJ   <- entropy_int(joint)
      mi    <- hI + hJ - hIJ
      if (abs(mi) < sqrt(.Machine$double.eps)) mi <- 0

      if (normalize) {
        emi <- .PairExpectedMI(tblI, tblJ)
        nmi <- .Rezero(mi / hBest, emi / hBest)
      } else {
        nmi <- mi / hBest
      }

      nmi_mat[pi, pj]   <- nmi_mat[pj, pi]   <- nmi
      hBest_mat[pi, pj] <- hBest_mat[pj, pi] <- hBest
    }
  }

  # ---- 6. Gaussian-weighted aggregation (-d analogue via sigma) ----
  maxHalf <- ceiling(truncate * max(sigma))
  scores  <- matrix(NA_real_, L, length(sigma),
                    dimnames = list(seq_len(L), sigma))

  for (i in seq_len(L)) {
    pi <- idx[i]
    if (patH[pi] == 0) next

    jLo <- max(1L, i - maxHalf)
    jHi <- min(L,  i + maxHalf)
    js  <- c(seq_len(i - jLo) + jLo - 1L,    # jLo..(i-1)
             seq_len(jHi - i) + i)             # (i+1)..jHi
    if (length(js) == 0L) next

    pjs     <- idx[js]
    nmiVals <- nmi_mat[cbind(rep(pi, length(pjs)), pjs)]
    hbVals  <- hBest_mat[cbind(rep(pi, length(pjs)), pjs)]

    valid <- !is.na(nmiVals) & hbVals > 0
    if (!any(valid)) next

    dists <- abs(js[valid] - i)
    nmiV  <- nmiVals[valid]
    hbV   <- hbVals[valid]

    for (s in seq_along(sigma)) {
      w  <- dnorm(dists, 0, sigma[s]) * hbV
      sw <- sum(w)
      scores[i, s] <- if (sw > 0) sum(w * nmiV) / sw else NA_real_
    }
  }

  # ---- 7. Block scan (-b / -p analogues) ----
  trim <- NULL
  if (!is.null(block_size)) {
    block_size <- as.integer(block_size)
    ref_scores <- scores[, 1L]    # smallest sigma: most local signal

    # Binary "noisy" label for each site
    is_noisy   <- !is.na(ref_scores) & ref_scores <= noise_level
    is_scored  <- !is.na(ref_scores)

    # Rolling sums over block_size using a rectangular kernel
    # filter() with sides = 2 centres the window; NAs at edges become NA
    kern <- rep(1, block_size)
    noisy_sum  <- suppressWarnings(
      as.numeric(stats::filter(as.numeric(is_noisy),  kern, sides = 2))
    )
    scored_sum <- suppressWarnings(
      as.numeric(stats::filter(as.numeric(is_scored), kern, sides = 2))
    )

    prop_noisy   <- ifelse(scored_sum > 0L, noisy_sum / scored_sum, NA_real_)
    in_bad_block <- !is.na(prop_noisy) & prop_noisy >= threshold

    # Trim from first noisy to last noisy site within each bad run,
    # matching PhyIN's "first to last conflicted" behaviour.
    trim      <- logical(L)
    run_open  <- FALSE
    run_start <- 0L
    last_noisy <- 0L

    for (i in seq_len(L)) {
      if (in_bad_block[i]) {
        if (!run_open) {
          run_open  <- TRUE
          run_start <- i
        }
        if (is_noisy[i]) last_noisy <- i
      } else {
        if (run_open) {
          if (last_noisy >= run_start) trim[run_start:last_noisy] <- TRUE
          run_open   <- FALSE
          last_noisy <- 0L
        }
      }
    }
    if (run_open && last_noisy >= run_start) {
      trim[run_start:last_noisy] <- TRUE
    }
  }

  structure(scores,
            pairNMI = nmi_mat,
            trim    = trim,
            class   = "LocalConcordance")
}


# ---- Internal: generalized expected MI for two multi-state characters ----
# Implements Vinh et al. (2010) Theorem 1 using the marginal hypergeometric
# distributions of each cell in the contingency table.
#' @importFrom stats dhyper
.PairExpectedMI <- function(a, b) {
  a <- a[a > 0L]
  b <- b[b > 0L]
  if (length(a) < 2L || length(b) < 2L) return(0)
  N     <- sum(a)
  log2N <- log2(N)
  emi   <- 0

  for (ai in a) {
    log2ai <- log2(ai)
    for (bj in b) {
      log2bj <- log2(bj)
      kmin <- max(1L, ai + bj - N)
      kmax <- min(ai, bj)
      if (kmin > kmax) next

      k  <- kmin:kmax
      pk <- dhyper(k, ai, N - ai, bj)
      emi <- emi + sum(pk * (k / N) * (log2(k) + log2N - log2ai - log2bj))
    }
  }
  emi
}


#' Plot local concordance
#'
#' Displays a multi-scale local concordance summary alongside (optionally) the
#' alignment itself, using base graphics.  Sites recommended for trimming
#' (the `"trim"` attribute, if present) are shaded in grey.
#'
#' @param lc A `LocalConcordance` matrix returned by [LocalConcordance()].
#' @param dataset Optional `phyDat` object; if provided, the alignment is drawn
#'   as a bottom panel.
#' @param cols Colour ramp (vector of hex codes) for concordance scores.
#'   Defaults to a red--white--blue diverging palette.
#' @param zlim Numeric vector of length 2 giving the score range mapped to the
#'   colour scale.  Defaults to the range of finite values in `lc`.
#' @param trim_col Colour used to shade trimmed regions (default light grey).
#'   Set to `NA` to suppress trim shading even when a `"trim"` attribute is
#'   present.
#' @param mar Numeric vector of length 4; base-graphics margins for each
#'   concordance panel.
#' @param \dots Additional arguments passed to [image()].
#'
#' @return Invisibly returns `lc`.
#'
#' @examples
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' dataset <- congreveLamsdellMatrices[[1]]
#' lc <- LocalConcordance(dataset, sigma = c(1, 1.5, 2, 5),
#'                        internal_gaps = FALSE,
#'                        block_size = 5L, noise_level = 0.1)
#' PlotLocalConcordance(lc, dataset)
#'
#' @seealso [LocalConcordance()] to compute scores.
#' @template MRS
#' @importFrom graphics axis image layout mtext par rect
#' @importFrom grDevices colorRampPalette
#' @importFrom TreeTools PhyDatToMatrix
#' @export
PlotLocalConcordance <- function(
  lc,
  dataset   = NULL,
  cols      = NULL,
  zlim      = NULL,
  trim_col  = "#CCCCCC60",
  mar       = c(0.5, 4, 0.5, 1),
  ...
) {
  if (is.null(cols)) {
    cols <- colorRampPalette(c("#D73027", "#FFFFFF", "#4575B4"))(256)
  }
  nCols  <- length(cols)
  nSigma <- ncol(lc)
  L      <- nrow(lc)
  sigma  <- as.numeric(colnames(lc))
  trim   <- attr(lc, "trim")

  if (is.null(zlim)) {
    zlim <- range(lc, na.rm = TRUE, finite = TRUE)
  }

  showAlignment <- !is.null(dataset)
  nPanels <- nSigma + showAlignment
  heights <- c(rep(1, nSigma), if (showAlignment) 3)
  layout(matrix(seq_len(nPanels)), heights = heights)

  # Helper: shade trim regions if present
  .shadeTrim <- function(trim, trim_col) {
    if (is.null(trim) || is.na(trim_col)) return(invisible(NULL))
    runs  <- rle(trim)
    ends  <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    for (r in which(runs$values)) {
      rect(starts[r] - 0.5, par("usr")[3],
           ends[r]   + 0.5, par("usr")[4],
           col = trim_col, border = NA)
    }
  }

  # -- Concordance panels (one per sigma, largest sigma at top) --
  for (s in rev(seq_len(nSigma))) {
    oPar <- par(mar = mar)
    on.exit(par(oPar), add = TRUE)

    vals <- lc[, s]
    scaled   <- (vals - zlim[1]) / diff(zlim)
    scaled[scaled < 0] <- 0
    scaled[scaled > 1] <- 1
    colIdx   <- pmin(nCols, floor(scaled * (nCols - 1)) + 1L)
    colIdx[is.na(colIdx)] <- NA_integer_

    plot.new()
    plot.window(xlim = c(0.5, L + 0.5), ylim = c(0, 1))
    rect(seq_len(L) - 0.5, 0, seq_len(L) + 0.5, 1,
         col = cols[colIdx], border = NA)
    .shadeTrim(trim, trim_col)
    mtext(bquote(sigma == .(sigma[s])), side = 2, las = 1, line = 1,
          cex = 0.8)
  }

  # -- Alignment panel --
  if (showAlignment) {
    par(mar = c(2.5, 4, 0, 1))
    charMat <- `mode<-`(PhyDatToMatrix(dataset), "numeric")
    nTaxa  <- nrow(charMat)
    nSites <- ncol(charMat)
    image(seq_len(nSites), seq_len(nTaxa), t(charMat),
          axes = FALSE, xlab = "", ylab = "")
    .shadeTrim(trim, trim_col)
    axis(1)
    mtext("Column", side = 1, line = 1.5, cex = 0.8)
    mtext("Taxon", side = 2, line = 1, cex = 0.8)
  }

  invisible(lc)
}
