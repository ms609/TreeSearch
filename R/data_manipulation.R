# Feasibility thresholds for MaddisonSlatkin exact computation.
# The split_count is the coefficient of x^floor(n/2) in the generating
# polynomial prod_i (1 + x + ... + x^{a_i}), capturing partition shape.
# Calibrated from worst-case (balanced) partition timing experiments
# using bitmask encoding (states at positions 2^(i-1)):
#   k=3: n=27 (9,9,9)   sc=75  0.97s safe;  n=31 (11,10,10) sc=96 1.32s marginal
#   k=4: n=13 (4,3,3,3) sc=50  0.36s safe;  n=15 (4,4,4,3)  sc=70 0.94s marginal
#   k=5: n=9  (2,2,2,2,1) sc=35 0.22s safe; n=10 (2,2,2,2,2) sc=51 0.49s
.MS_SC_THRESHOLD <- c(Inf, Inf, 75L, 50L, 35L)

.MSSplitCount <- function(state_counts) {
  counts <- state_counts[state_counts > 0L]
  if (!length(counts)) return(0L)
  n <- sum(counts)
  if (n <= 2L) return(1L)
  target <- n %/% 2L
  poly <- 1.0
  for (ci in counts) {
    new_len  <- min(length(poly) + ci, target + 1L)
    new_poly <- numeric(new_len)
    for (j in seq_len(new_len)) {
      lo <- max(1L, j - ci)
      hi <- min(j,  length(poly))
      if (lo <= hi) new_poly[j] <- sum(poly[lo:hi])
    }
    poly <- new_poly
  }
  if (target + 1L <= length(poly)) poly[target + 1L] else 0.0
}

#' Prepare data for Profile Parsimony
#' 
#' Calculates profiles for each character in a dataset.
#' Characters with 2 informative states (i.e. states present in more than one
#' taxon) use the exact formula of Carter _et al._ (1990).
#' Characters with 3 or more informative states use the recursive algorithm of
#' Maddison & Slatkin (1991), falling back to a Monte Carlo approximation for
#' large or complex characters.
#' 
#' Characters are simplified where necessary, with a warning:
#' - inapplicable tokens will be replaced with the ambiguous token
#'    (i.e. `-` \ifelse{html}{\out{&rarr;}}{\eqn{\rightarrow}{-->}} `?`);
#' - Ambiguous tokens will be treated as fully ambiguous
#'   (i.e. `{02}` \ifelse{html}{\out{&rarr;}}{\eqn{\rightarrow}{-->}} `?`)
#' 
#' @param dataset dataset of class \code{phyDat}
#' @param approx Character string controlling how profile information amounts
#'   are computed for multi-state characters with many tips.
#'   `"auto"` (default) uses the exact Maddison & Slatkin calculation when
#'   feasible, falling back to a Monte Carlo approximation for large or
#'   complex characters.
#'   `"mc"` always uses the Monte Carlo approximation;
#'   `"exact"` always uses the exact calculation (may be very slow).
#' @param n_mc Integer; number of Monte Carlo samples for the MC
#'   approximation.  Default 100 000.
#'
#' @return An object of class `phyDat`, with additional attributes.
#' `PrepareDataProfile` adds the attributes:
#' 
#'   - `info.amounts`: details the information represented by each 
#'     character when subject to N additional steps.
#'   
#'   - `informative`: logical specifying which characters contain any
#'     phylogenetic information.
#'   
#'   - `bootstrap`: The character vector 
#'     \code{c("info.amounts", "split.sizes")}, indicating attributes to sample
#'      when bootstrapping the dataset (e.g. in Ratchet searches).
#'
#' `PrepareDataIW` adds the attribute:
#' 
#'  - \code{min.length}: The minimum number of steps that must be present in each
#'    transformation series.
#' @examples 
#' data("congreveLamsdellMatrices")
#' dataset <- congreveLamsdellMatrices[[42]]
#' PrepareDataProfile(dataset)
#' @author Martin R. Smith; written with reference to 
#' `phangorn:::prepareDataFitch()`
#' @importFrom cli cli_alert cli_alert_warning
#' @importFrom fastmatch %fin%
#' @family profile parsimony functions
#' @encoding UTF-8
#' @export
PrepareDataProfile <- function (dataset, approx = "auto", n_mc = 100000L) {
  if ("info.amounts" %fin% names(attributes(dataset))) {
    # Already prepared
    return(dataset)
  }
  at <- attributes(dataset)
  nLevel <- length(at[["levels"]])
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  index <- at[["index"]]
  allLevels <- as.character(at[["allLevels"]])
  
  contSums <- rowSums(cont)
  qmLevel <- which(contSums == ncol(cont))
  
  if (length(qmLevel) == 0) {
    attr(dataset, "contrast") <- rbind(attr(dataset, "contrast"), 1)
    attr(dataset, "allLevels") <- c(attr(dataset, "allLevels"), "{?}")
    qmLevel <- length(allLevels) + 1L
  }
  
  ambigs <- which(contSums > 1L & contSums < ncol(cont))
  inappLevel <- which(colnames(cont) == "-")
  if (length(inappLevel) != 0L) {
    cli_alert("Inapplicable tokens treated as ambiguous for profile parsimony")
    inappLevel <- which(apply(unname(cont), 1, identical,
                              as.double(colnames(cont) == "-")))
    dataset[] <- lapply(dataset, function (i) {
      i[i %fin% inappLevel] <- qmLevel
      i
    })
  }
  
  if (length(ambigs) != 0L) {
    dataset[] <- lapply(dataset, function (i) {
        i[i %fin% ambigs] <- qmLevel
        i
      })
  }
  
  # Build pattern matrix: rows = patterns (unique characters), cols = tips
  nPattern <- max(index)
  mataset <- matrix(unlist(dataset, recursive = FALSE, use.names = FALSE),
                    nPattern)
  # Transpose to: rows = tips, cols = patterns (matching .RemoveExtraTokens)
  mataset <- t(mataset)
  
  # --- Strip singletons ---
  maxInformative <- 0L
  
  for (j in seq_len(ncol(mataset))) {
    col <- mataset[, j]
    nonAmbig <- col[col != qmLevel[1]]
    if (length(nonAmbig) == 0L) next
    
    tab <- table(nonAmbig)
    informative <- tab > 1L
    nInf <- sum(informative)
    
    # Convert singletons to ambiguous
    singletonTokens <- as.integer(names(tab[!informative]))
    if (length(singletonTokens) > 0L) {
      mataset[mataset[, j] %in% singletonTokens, j] <- qmLevel[1]
    }
    
    maxInformative <- max(maxInformative, nInf)
  }

  
  if (maxInformative < 2L) {
    cli_alert("No informative characters in `dataset`.")
    # Construct empty phyDat manually (avoids [.phyDat issues with 0 columns)
    dataset[] <- lapply(dataset, function(x) integer(0))
    attr(dataset, "info.amounts") <- double(0)
    attr(dataset, "weight") <- integer(0)
    attr(dataset, "nr") <- 0L
    attr(dataset, "index") <- integer(0)
    return(dataset)
  }
  
  # --- Recompress: normalize tokens to 1..k, AMBIG ---
  AMBIG_TOKEN <- maxInformative + 1L
  
  for (j in seq_len(ncol(mataset))) {
    col <- mataset[, j]
    nonAmbig <- sort(unique(col[col != qmLevel[1]]))
    
    newCol <- rep(AMBIG_TOKEN, length(col))
    for (i in seq_along(nonAmbig)) {
      newCol[col == nonAmbig[i]] <- i
    }
    mataset[, j] <- newCol
  }
  
  # --- Deduplicate patterns ---
  dupCols <- duplicated(t(mataset))
  kept <- which(!dupCols)
  copies <- lapply(kept, function (i) {
    i + which(apply(mataset[, -seq_len(i), drop = FALSE], 2, identical, mataset[, i]))
  })
  firstOccurrence <- seq_len(dim(mataset)[2])
  for (i in seq_along(copies)) {
    firstOccurrence[copies[[i]]] <- kept[i]
  }
  
  cipher <- seq_len(max(kept))
  cipher[kept] <- order(kept)
  index <- cipher[firstOccurrence][index]
  
  mataset <- mataset[, !dupCols, drop = FALSE]
  dataset[] <- lapply(seq_len(length(dataset)), function (i) mataset[i, ])
  
  # --- Compute StepInformation per unique pattern ---
  info <- lapply(seq_along(mataset[1, ]), function (i) 
    StepInformation(mataset[, i], ambiguousTokens = AMBIG_TOKEN,
                    approx = approx, n_mc = n_mc))
  
  
  maxSteps <- max(vapply(info,
                         function (i) max(as.integer(names(i))),
                         integer(1)))
  info <- vapply(info,
                 function (x) {
                    ret <- setNames(double(maxSteps), seq_len(maxSteps))
                    x <- x[setdiff(names(x), "0")]
                    if (length(x)) {
                      ret[names(x)] <- max(x) - x
                    }
                    ret
                  }, double(maxSteps))
  if (is.null(dim(info))) {
    dim(info) <- c(1L, length(info))
  }
  attr(dataset, "index") <- index
  weight <- as.integer(table(index))
  attr(dataset, "weight") <- weight
  attr(dataset, "nr") <- length(weight)
  attr(dataset, "info.amounts") <- info
  attr(dataset, "informative") <- colSums(info) > 0
  
  # Dynamic contrast matrix: k states + ambiguous
  k <- maxInformative
  lvls <- as.character(seq_len(k))
  contMatrix <- rbind(diag(k), rep(1L, k))
  dimnames(contMatrix) <- list(NULL, lvls)
  
  attr(dataset, "levels") <- lvls
  attr(dataset, "allLevels") <- c(lvls, "?")
  attr(dataset, "contrast") <- contMatrix
  attr(dataset, "nc") <- as.integer(k)
  
  if (!any(attr(dataset, "bootstrap") == "info.amounts")) {
    attr(dataset, "bootstrap") <- c(attr(dataset, "bootstrap"), "info.amounts")
  }
  
  dataset
}


#' @describeIn PrepareDataProfile Prepare data for implied weighting
#' @export
PrepareDataIW <- function (dataset) {
  attr(dataset, "min.length") <- MinimumLength(dataset, compress = TRUE)
  
  # Return:
  dataset
}
