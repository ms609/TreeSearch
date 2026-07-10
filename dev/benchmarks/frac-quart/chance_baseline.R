# Validate the exact hypergeometric chance baseline used by
# QuartetConcordance(unit = "trit", normalize = TRUE) against a Monte-Carlo
# tip-shuffle oracle, at the level of a single (split, character, state-pair).
#
# Under the fixed-marginal null the character's tokens are reassigned across the
# t scored leaves at random, holding the pair's state counts (n_i, n_j) and the
# split's scored side-A size M fixed.  The cell vector (p,q,r,s) is then
# trivariate-hypergeometric.  We need the EXPECTED per-pair pool contributions
#   m       = min(w_c, w_k)          (shared-information pooling weight)
#   m*A/w_k (edge-numerator contribution)
#   m*A/w_c (char-numerator contribution)
# because w_k = (mA-1)+(tP-mA-1)+ -- and hence m -- are THEMSELVES random for
# multistate pairs (mA = p + r varies).  Binary characters are the special case
# mA == M, tP == t, so w_k is fixed and only A varies (the clean
# ClusteringConcordance case).  This mirrors `.ExpectedTrit()` in R/Concordance.R.

pos <- function(z) { z[z < 0] <- 0; z }

# Observed per-pair triple for one realised cell vector.
pairTriple <- function(p, q, r, s) {
  nI <- p + q; nJ <- r + s
  mA <- p + r; tP <- nI + nJ
  A  <- pos(p - 1) * pos(s - 1) + pos(q - 1) * pos(r - 1)
  wc <- pos(nI - 1) * pos(nJ - 1)
  wk <- pos(mA - 1) * pos(tP - mA - 1)
  m  <- min(wc, wk)
  c(m    = m,
    mAwk = if (wk > 0) m * A / wk else 0,
    mAwc = if (wc > 0) m * A / wc else 0)
}

# Exact expectation over the trivariate hypergeometric (the R kernel's method).
expectedTriple <- function(nI, nJ, M, t) {
  nOther <- t - nI - nJ
  stopifnot(nOther >= 0, M >= 0, M <= t)
  logC <- function(n, k) if (k < 0 || k > n) -Inf else lchoose(n, k)
  denom <- lchoose(t, M)
  acc <- c(m = 0, mAwk = 0, mAwc = 0)
  for (p in 0:min(nI, M)) {
    for (r in 0:min(nJ, M - p)) {
      oA <- M - p - r
      if (oA < 0 || oA > nOther) next
      prob <- exp(logC(nI, p) + logC(nJ, r) + logC(nOther, oA) - denom)
      if (prob <= 0) next
      acc <- acc + prob * pairTriple(p, nI - p, r, nJ - r)
    }
  }
  acc
}

# Monte-Carlo oracle: shuffle state labels across t leaves, M fixed on side A.
mcTriple <- function(nI, nJ, M, t, nShuf = 2e5) {
  sideA <- c(rep(TRUE, M), rep(FALSE, t - M))
  acc <- c(m = 0, mAwk = 0, mAwc = 0)
  for (i in seq_len(nShuf)) {
    lab <- sample(rep(1:3, c(nI, nJ, t - nI - nJ)))   # 1 = i, 2 = j, 3 = other
    p <- sum(lab == 1 & sideA); q <- sum(lab == 1 & !sideA)
    r <- sum(lab == 2 & sideA); s <- sum(lab == 2 & !sideA)
    acc <- acc + pairTriple(p, q, r, s)
  }
  acc / nShuf
}

set.seed(1)
cases <- list(
  c(nI = 4, nJ = 4, M = 4, t = 8),   # binary-like (nOther = 0): w_k fixed
  c(nI = 3, nJ = 5, M = 4, t = 8),
  c(nI = 3, nJ = 3, M = 4, t = 9),   # multistate (nOther > 0): w_k random
  c(nI = 4, nJ = 2, M = 5, t = 10),
  c(nI = 5, nJ = 4, M = 6, t = 12),
  c(nI = 2, nJ = 2, M = 3, t = 7),
  c(nI = 6, nJ = 3, M = 5, t = 15)
)

worst <- 0
for (cs in cases) {
  ex <- expectedTriple(cs["nI"], cs["nJ"], cs["M"], cs["t"])
  mc <- mcTriple(cs["nI"], cs["nJ"], cs["M"], cs["t"])
  worst <- max(worst, max(abs(ex - mc)))
  cat(sprintf("nI=%d nJ=%d M=%2d t=%2d  exact=(%.4f,%.4f,%.4f)  mc=(%.4f,%.4f,%.4f)\n",
              cs["nI"], cs["nJ"], cs["M"], cs["t"],
              ex[1], ex[2], ex[3], mc[1], mc[2], mc[3]))
}
cat(sprintf("\nWorst exact-vs-MC discrepancy: %.5f\n", worst))
stopifnot(worst < 0.05)
cat("PASS: exact hypergeometric baseline matches the MC tip-shuffle oracle.\n")
