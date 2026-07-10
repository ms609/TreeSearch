# Parity guard for the C++ chance-baseline kernels `quartet_expect` and
# `trit_expect` (src/concordance_expect.cpp), which compute the exact expected
# concordant/decisive quartet counts (raw) and the expected trit pools
# (E[m], E[m*A/wk], E[m*A/wc]) under the fixed-marginal (hypergeometric) null.
#
# The kernels replaced an earlier pure-R implementation; this script keeps a
# self-contained R reference (the same trivariate-hypergeometric double sum) and
# checks the C++ output matches it to machine precision across random binary,
# multistate and missing-data matrices.  Run after any change to the kernels.

suppressMessages(pkgload::load_all(
  "C:/Users/pjjg18/GitHub/worktrees/TreeSearch/frac-quart-cpp",
  quiet = TRUE, recompile = TRUE))
suppressMessages(library(TreeTools))

# Build logiSplits + charInt exactly as QuartetConcordance() does internally.
build_inputs <- function(tree, dat) {
  tl <- intersect(TipLabels(tree), names(dat))
  dat <- dat[tl, drop = FALSE]
  splits <- as.Splits(tree, dat)
  logiSplits <- vapply(seq_along(splits), function(i) as.logical(splits[[i]]),
                       logical(NTip(dat)))
  contrast <- attr(dat, "contrast"); lv <- attr(dat, "allLevels")
  isInapp <- lv == "-"
  isAmbig <- rowSums(contrast[, colnames(contrast) != "-", drop = FALSE]) > 1
  isGrouping <- !isAmbig & !isInapp
  gCols <- apply(contrast[isGrouping, , drop = FALSE] > 0, 1, which)
  l2i <- rep(NA_integer_, length(lv)); l2i[isGrouping] <- as.integer(gCols)
  ch <- PhyDatToMatrix(dat)
  charInt <- array(l2i[match(ch, lv)], dim = dim(ch), dimnames = dimnames(ch))
  list(logiSplits = logiSplits, charInt = charInt)
}

# ---- self-contained R reference (trivariate-hypergeometric double sum) ----
pos1 <- function(z) if (z < 0) 0 else z
choose2 <- function(z) if (z < 2) 0 else z * (z - 1) / 2

ref_pair <- function(nI, nJ, M, t) {           # returns c(Econc, Edec, Em, Ewk, Ewc)
  nOther <- t - nI - nJ; tP <- nI + nJ
  wc <- pos1(nI - 1) * pos1(nJ - 1)
  logDen <- lchoose(t, M)
  acc <- c(0, 0, 0, 0, 0)
  for (p in 0:min(nI, M)) for (r in 0:min(nJ, M - p)) {
    oA <- M - p - r
    if (oA < 0 || oA > nOther) next
    prob <- exp(lchoose(nI, p) + lchoose(nJ, r) + lchoose(nOther, oA) - logDen)
    if (prob <= 0) next
    q <- nI - p; s <- nJ - r; mA <- p + r
    conc <- choose2(p) * choose2(s) + choose2(q) * choose2(r)
    dec <- conc + p * q * r * s
    A <- pos1(p - 1) * pos1(s - 1) + pos1(q - 1) * pos1(r - 1)
    wk <- pos1(mA - 1) * pos1(tP - mA - 1); m <- min(wc, wk)
    acc <- acc + prob * c(conc, dec, m, if (wk > 0) m * A / wk else 0,
                          if (wc > 0) m * A / wc else 0)
  }
  acc
}

ref_expect <- function(charInt, logiSplits) {
  nSplit <- ncol(logiSplits); nChar <- ncol(charInt)
  out <- list(eConc = matrix(0, nSplit, nChar), eDec = matrix(0, nSplit, nChar),
              eM = matrix(0, nSplit, nChar), eWk = matrix(0, nSplit, nChar),
              eWc = matrix(0, nSplit, nChar))
  for (ci in seq_len(nChar)) {
    col <- charInt[, ci]; scored <- !is.na(col)
    states <- sort(unique(col[scored])); nStates <- length(states)
    if (nStates < 2L) next
    tc <- sum(scored); mSideA <- colSums(logiSplits & scored)
    uM <- unique(mSideA); idx <- match(mSideA, uM)
    cnt <- tabulate(match(col[scored], states), nStates)
    for (a in seq_len(nStates - 1L)) for (b in seq(a + 1L, nStates)) {
      e <- vapply(uM, function(M) ref_pair(cnt[a], cnt[b], M, tc), double(5))
      out$eConc[, ci] <- out$eConc[, ci] + e[1, idx]
      out$eDec[, ci]  <- out$eDec[, ci]  + e[2, idx]
      out$eM[, ci]    <- out$eM[, ci]    + e[3, idx]
      out$eWk[, ci]   <- out$eWk[, ci]   + e[4, idx]
      out$eWc[, ci]   <- out$eWc[, ci]   + e[5, idx]
    }
  }
  out
}

set.seed(7)
worst <- 0
for (rep in 1:10) {
  nTip <- sample(9:44, 1)
  tree <- RandomTree(nTip, root = FALSE)
  m <- replicate(sample(5:30, 1), {
    ns <- sample(2:4, 1)
    x <- as.character(sample(0:(ns - 1), nTip, replace = TRUE))
    if (runif(1) < 0.3) x[sample(nTip, sample(1:3, 1))] <- "?"
    x
  })
  rownames(m) <- paste0("t", seq_len(nTip))
  io <- build_inputs(tree, MatrixToPhyDat(m))
  ls <- io$logiSplits; ci <- io$charInt
  R <- ref_expect(ci, ls)
  cq <- quartet_expect(ls, ci); ct <- trit_expect(ls, ci)
  d <- max(abs(cq$concordant - R$eConc), abs(cq$decisive - R$eDec),
           abs(ct$numEdge - R$eWk), abs(ct$numChar - R$eWc),
           abs(ct$denM - R$eM))
  worst <- max(worst, d)
  cat(sprintf("rep %2d: nTip=%2d nChar=%2d  max|C++ - R| = %.2e\n",
              rep, nTip, ncol(ci), d))
}
cat(sprintf("\nWorst |C++ - R| over all kernels: %.3e\n", worst))
stopifnot(worst < 1e-9)
cat("PASS: C++ expectation kernels match the R reference to machine precision.\n")
