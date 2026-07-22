# Fractional per-character weights.
#
# TreeSearch's C++ scoring engine stores per-pattern weights as `int`.
# Without intervention, a phyDat with `attr(dat, "weight") <- c(0.5, 1.7)`
# would silently truncate to `c(0L, 1L)` at the Rcpp boundary, dropping
# 50% of the first character's contribution and 41% of the second's.
#
# `.ScaleWeight()` converts a fractional weight vector to integer with a
# documented scale factor (default 2*2*3*3*5*7 = 1260, ~0.001 precision).
# Integerweights pass through unchanged so the function is a no-op for the
# common case.
#
# The TreeLength value returned by the scoring engine is then in units of
# (steps * scale), so users comparing across runs with fractional weights
# should divide by `getOption("TreeSearch.fractional.scale", 1260L)`
# (or rely on within-run ranking, which is unaffected).

#' @keywords internal
.ScaleWeight <- function(weight) {
  if (length(weight) == 0L) {
    # Return:
    return(integer(0L))
  }
  # Reject values that would corrupt the integer weight passed to C++: a
  # negative weight reaches the scorer as a negative `int` (undefined
  # behaviour), and NA/NaN/Inf otherwise surface only as an opaque
  # "missing value where TRUE/FALSE needed" from the overflow guard below.
  if (any(!is.finite(weight)) || any(weight < 0)) {
    stop("`weight` must contain only finite, non-negative values.",
         call. = FALSE)
  }
  if (is.integer(weight)) {
    # Return:
    weight
  } else if (all(weight == as.integer(weight))) {
    # Already-integer values stored as double: cast and return without scaling.
    # Return:
    as.integer(weight)
  } else {
    scale <- as.integer(getOption("TreeSearch.fractional.scale", 1260L))
    if (scale < 1L) scale <- 1L
    scaled <- as.integer(round(weight * scale))
    # Guard against under-rounded weights becoming zero: a weight of zero
    # would silently drop the character. Floor at 1 unless the user
    # genuinely supplied a zero weight (preserved as 0L).
    keep <- weight > 0
    scaled[keep & scaled < 1L] <- 1L
    # Guard: the C++ resampling routines expand each pattern `weight[p]`
    # times into a flat index vector whose length is cast to `int`.  If
    # sum(weights) > .Machine$integer.max the cast overflows to a negative
    # value and the subsequent array access is undefined behaviour (segfault).
    total <- sum(as.double(scaled))
    if (total > .Machine$integer.max) {
      stop(
        "Total scaled weight (",
        format(round(total), big.mark = ",", scientific = FALSE),
        ") exceeds .Machine$integer.max.\n",
        "Reduce options(\"TreeSearch.fractional.scale\") from ", scale,
        " to a smaller value, or set integer weights directly.",
        call. = FALSE
      )
    }
    # Return:
    scaled
  }
}
