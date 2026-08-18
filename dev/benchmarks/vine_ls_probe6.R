# ---------------------------------------------------------------------------
# Probe 6 -- the one actionable consequence.
#
# The gradient half of the Vine idea fails (probes 1-5).  The DECODER half
# does not: perturbing the input distance matrix and re-running NJ is a
# surprisingly effective, and very cheap, move operator -- it found the exact
# NNI+SPR optimum (RF = 0) on every target, 62x faster at n=40/low conflict.
#
# So: does a cheap perturbed-NJ pre-pass reduce the total cost of
# LeastSquaresTree() to the same optimum?
#
#   (A) current default        : LeastSquaresTree(D)                 [NJ start]
#   (B) pre-pass then search   : K perturbed-NJ decodes -> best -> search
#
# Cost of (B) includes the pre-pass.  Both must reach the same RSS to count.
# ---------------------------------------------------------------------------

suppressMessages({
  library(TreeSearch, lib.loc = if (dir.exists(Sys.getenv("TS_LIB", ".agent-p0"))) Sys.getenv("TS_LIB", ".agent-p0") else .libPaths())
  library(ape)
})
SCR <- Sys.getenv("SCRATCH", unset = ".")
METHOD <- "nnls"; K <- 150L; EPS <- 0.15

MakeTarget <- function(nTip, m, k = 12L, seed = 7L) {
  set.seed(seed)
  t0 <- TreeTools::RandomTree(nTip, root = FALSE)
  t0[["tip.label"]] <- sprintf("t%02d", seq_len(nTip))
  trees <- c(list(t0), lapply(seq_len(k - 1L), function(i) phangorn::rSPR(t0, moves = m)))
  labs <- sort(t0[["tip.label"]])
  P <- lapply(trees, function(tr) {
    tr[["edge.length"]] <- rep(1, nrow(tr[["edge"]])); cophenetic(tr)[labs, labs] })
  D <- Reduce(`+`, P) / length(P); diag(D) <- 0
  list(D = D, labs = labs)
}
RealTarget <- function(dsName) {
  mpts <- readRDS(file.path(SCR, paste0("mpt-", dsName, ".rds")))
  labs <- sort(mpts[[1]][["tip.label"]])
  P <- lapply(mpts, function(tr) {
    tr[["edge.length"]] <- rep(1, nrow(tr[["edge"]])); cophenetic(tr)[labs, labs] })
  D <- Reduce(`+`, P) / length(P); diag(D) <- 0
  list(D = D, labs = labs)
}

Compare <- function(label, tgt) {
  D <- tgt$D; nTip <- nrow(D); pairIdx <- upper.tri(D); sdD <- sd(D[pairIdx])
  Decode <- function(M) { tr <- ape::nj(as.dist(M))
    tr[["edge.length"]][tr[["edge.length"]] < 0] <- 0; tr }
  Perturb <- function() { z <- matrix(0, nTip, nTip); z[pairIdx] <- rnorm(sum(pairIdx))
    Mp <- D + EPS * sdD * (z + t(z)); Mp[Mp < 0] <- 0; diag(Mp) <- 0
    dimnames(Mp) <- dimnames(D); Mp }
  Refit <- function(tr) { f <- try(LeastSquaresFit(tr, D, method = METHOD), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(attr(f, "RSS")) }

  tA <- system.time(A <- LeastSquaresTree(D, method = METHOD))[["elapsed"]]
  rA <- as.numeric(attr(A, "RSS"))

  tB <- system.time({
    bv <- Inf; bt <- NULL
    for (i in seq_len(K)) { tr <- Decode(Perturb()); v <- Refit(tr)
      if (is.finite(v) && v < bv) { bv <- v; bt <- tr } }
    seed0 <- if (is.null(bt)) NULL else bt
    B <- LeastSquaresTree(D, tree = seed0, method = METHOD)
  })[["elapsed"]]
  rB <- as.numeric(attr(B, "RSS"))

  cat(sprintf("%-26s n=%2d | (A) default %9.4f in %6.2f s | pre-pass best %9.4f | (B) total %9.4f in %6.2f s | %s %.1fx\n",
              label, nTip, rA, tA, bv, rB, tB,
              if (rB < rA - 1e-9) "BETTER," else if (rB > rA + 1e-9) "WORSE," else "same,",
              tA / tB))
  data.frame(label = label, nTip = nTip, defaultRSS = rA, defaultSec = tA,
             prepassBest = bv, seededRSS = rB, seededSec = tB,
             speedup = tA / tB)
}

set.seed(303)
res <- rbind(
  Compare("sim n=20 conflict m=2", MakeTarget(20, 2)),
  Compare("sim n=20 conflict m=8", MakeTarget(20, 8)),
  Compare("sim n=40 conflict m=2", MakeTarget(40, 2)),
  Compare("sim n=40 conflict m=8", MakeTarget(40, 8)),
  Compare("sim n=60 conflict m=4", MakeTarget(60, 4)),
  Compare("real Longrich2010",     RealTarget("Longrich2010")))
saveRDS(res, file.path(SCR, "vine-ls-probe6.rds"))
cat("\n"); print(res, row.names = FALSE, digits = 5)
