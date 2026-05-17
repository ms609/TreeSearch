# Fractional per-character weights.
#
# TreeSearch's C++ scoring engine stores per-pattern weights as `int`.
# Without intervention, a phyDat with `attr(dat, "weight") <- c(0.5, 1.7)`
# would silently truncate to `c(0L, 1L)` at the Rcpp boundary, dropping
# 50% of the first character's contribution and 41% of the second's.
#
# `.ScaleWeight()` converts a fractional weight vector to integer with a
# documented scale factor (default 1000, i.e. 0.001 precision). Integer
# weights pass through unchanged so the function is a no-op for the
# common case.
#
# The TreeLength value returned by the scoring engine is then in units of
# (steps * scale), so users comparing across runs with fractional weights
# should divide by `getOption("TreeSearch.fractional.scale", 1000L)`
# (or rely on within-run ranking, which is unaffected).

#' @keywords internal
.ScaleWeight <- function(weight) {
  if (length(weight) == 0L) {
    # Return:
    integer(0L)
  } else if (is.integer(weight)) {
    # Return:
    weight
  } else if (all(weight == as.integer(weight))) {
    # Already-integer values stored as double: cast and return without scaling.
    # Return:
    as.integer(weight)
  } else {
    scale <- as.integer(getOption("TreeSearch.fractional.scale", 1000L))
    if (scale < 1L) scale <- 1L
    scaled <- as.integer(round(weight * scale))
    # Guard against under-rounded weights becoming zero: a weight of zero
    # would silently drop the character. Floor at 1 unless the user
    # genuinely supplied a zero weight (preserved as 0L).
    keep <- weight > 0
    scaled[keep & scaled < 1L] <- 1L
    # Return:
    scaled
  }
}
