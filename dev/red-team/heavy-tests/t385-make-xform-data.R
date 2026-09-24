# Shared synthetic XFORM dataset generator for the T-385 scripts.
#
# nBlock hierarchy blocks, each a controlling primary + nSec secondaries, plus
# some free (non-hierarchy) characters.  The gain/loss asymmetry that makes XFORM
# rooting-sensitive needs nSec >= 1; the free characters carry no gaps, so
# `has_na` is driven only by the hierarchy blocks (the has_na = FALSE branch of
# the T-374b analysis, i.e. the Sankoff term stays out of the convergence sweep).

MakeXformData <- function(nTip = 36L, nBlock = 6L, nSec = 2L, nFree = 10L,
                          seed = 1L) {
  set.seed(seed)
  tipNames <- paste0("t", seq_len(nTip))
  cols <- list()
  hierArgs <- list()

  for (b in seq_len(nBlock)) {
    primary <- sample(c("0", "1"), nTip, replace = TRUE)
    priIdx <- length(cols) + 1L
    cols[[length(cols) + 1L]] <- primary
    secIdx <- integer(nSec)
    for (s in seq_len(nSec)) {
      # Secondary is inapplicable exactly where the primary is absent.
      sec <- ifelse(primary == "0", "-",
                    sample(c("0", "1"), nTip, replace = TRUE))
      cols[[length(cols) + 1L]] <- sec
      secIdx[s] <- length(cols)
    }
    hierArgs[[as.character(priIdx)]] <- secIdx
  }

  for (f in seq_len(nFree)) {
    cols[[length(cols) + 1L]] <- sample(c("0", "1"), nTip, replace = TRUE)
  }

  mat <- do.call(cbind, cols)
  dimnames(mat) <- list(tipNames, NULL)
  ds <- phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                         ambiguity = "?")
  list(dataset = ds, hierarchy = do.call(CharacterHierarchy, hierArgs),
       nSec = nSec, nBlock = nBlock)
}
