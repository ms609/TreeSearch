# ---------------------------------------------------------------------------
# Probe 7 -- turn the ORACLE timing into an honest, self-terminating procedure.
#
# Probe 5 reported "reached the optimum at draw 78, 0.89 s ==> 62x faster".
# That was measured KNOWING the answer.  A real sampler must decide when to
# stop.  Here: stop after `patience` consecutive draws with no improvement.
#
# Measured in CANDIDATE EVALUATIONS, not seconds -- both the sampler and
# NNI+SPR pay one NNLS refit per candidate, so evaluations are the honest
# common currency, and the number is immune to CPU contention.
#
# Two arms, because eps = 0.15 in probes 4-6 was itself chosen with hindsight
# from a sweep on these same targets:
#    fixed  : eps = 0.15                     (tuned)
#    random : eps ~ U(0.05, 0.30) per draw   (untuned, no oracle knowledge)
#
# One long chain per replicate; the stopping rule for every `patience` is then
# derived offline from the recorded trajectory.
# ---------------------------------------------------------------------------

suppressMessages({
  library(TreeSearch, lib.loc = if (dir.exists(Sys.getenv("TS_LIB", ".agent-p0"))) Sys.getenv("TS_LIB", ".agent-p0") else .libPaths())
  library(ape)
})
SCR <- Sys.getenv("SCRATCH", unset = ".")
METHOD <- "nnls"
PATIENCE <- c(50L, 200L, 800L, 2000L)

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

Run <- function(label, tgt, nRep = 10L, nDraw = 3000L) {
  D <- tgt$D; nTip <- nrow(D); pairIdx <- upper.tri(D); sdD <- sd(D[pairIdx])
  cat(sprintf("\n=========== %s (%d tips) ===========\n", label, nTip))
  Decode <- function(M) { tr <- ape::nj(as.dist(M))
    tr[["edge.length"]][tr[["edge.length"]] < 0] <- 0; tr }
  Perturb <- function(e) { z <- matrix(0, nTip, nTip); z[pairIdx] <- rnorm(sum(pairIdx))
    Mp <- D + e * sdD * (z + t(z)); Mp[Mp < 0] <- 0; diag(Mp) <- 0
    dimnames(Mp) <- dimnames(D); Mp }
  Refit <- function(tr) { f <- try(LeastSquaresFit(tr, D, method = METHOD), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(attr(f, "RSS")) }

  nj <- Decode(D); rssNJ <- Refit(nj)
  tS <- system.time(win <- LeastSquaresTree(D, method = METHOD))[["elapsed"]]
  rssS <- as.numeric(attr(win, "RSS"))
  tFit <- system.time(for (i in 1:100) Refit(nj))[["elapsed"]] / 100
  evalsSearch <- tS / tFit
  cat(sprintf("NJ %.4f | NNI+SPR %.4f | %.1f s = ~%.0f refit-equivalents (%.1f ms/refit)\n",
              rssNJ, rssS, tS, evalsSearch, 1000 * tFit))

  # ---- long chains, both eps arms ----
  chains <- list()
  for (arm in c("fixed", "random")) {
    for (r in seq_len(nRep)) {
      set.seed(1000L * match(arm, c("fixed", "random")) + r)
      v <- numeric(nDraw)
      for (i in seq_len(nDraw)) {
        e <- if (arm == "fixed") 0.15 else runif(1, 0.05, 0.30)
        x <- Refit(Decode(Perturb(e)))
        v[i] <- if (is.finite(x)) x else Inf
      }
      chains[[length(chains) + 1L]] <- list(arm = arm, rep = r, v = v)
    }
  }
  ref <- min(c(rssS, rssNJ, unlist(lapply(chains, function(c) min(c$v)))))
  cat(sprintf("reference best-known RSS = %.6f  (%s)\n", ref,
              if (ref < rssS - 1e-9) "sampler BEAT NNI+SPR" else "= NNI+SPR"))

  # ---- derive the stopping rule offline ----
  Stop <- function(v, patience) {
    best <- Inf; sinceImp <- 0L
    for (i in seq_along(v)) {
      if (v[i] < best - 1e-12) { best <- v[i]; sinceImp <- 0L } else sinceImp <- sinceImp + 1L
      if (sinceImp >= patience) return(c(evals = i, rss = best))
    }
    c(evals = length(v), rss = best)     # never triggered: censored
  }
  out <- do.call(rbind, lapply(chains, function(ch) {
    do.call(rbind, lapply(PATIENCE, function(p) {
      s <- Stop(ch$v, p)
      data.frame(arm = ch$arm, rep = ch$rep, patience = p,
                 evals = s[["evals"]], rss = s[["rss"]],
                 hit = s[["rss"]] <= ref + 1e-9,
                 censored = s[["evals"]] == length(ch$v))
    }))
  }))
  agg <- do.call(rbind, lapply(split(out, list(out$arm, out$patience), drop = TRUE),
    function(d) data.frame(arm = d$arm[1], patience = d$patience[1],
      pHit = mean(d$hit), medEvals = median(d$evals),
      medExcess = median(d$rss - ref), pCensored = mean(d$censored),
      evalRatio = evalsSearch / median(d$evals))))
  agg <- agg[order(agg$arm, agg$patience), ]
  cat(sprintf("\nNNI+SPR reference cost = ~%.0f evaluations\n", evalsSearch))
  print(agg, row.names = FALSE, digits = 4)
  cat("  pHit = P(returns the best-known tree); evalRatio = NNI+SPR evals / sampler evals\n")
  cat("  pCensored > 0 means the chain hit its 'nDraw' cap before stalling: evals is a LOWER bound\n")
  agg$label <- label; agg$nTip <- nTip; agg$evalsSearch <- evalsSearch
  agg$msPerFit <- 1000 * tFit
  agg
}

res <- rbind(
  Run("sim n=20 conflict m=8", MakeTarget(20, 8), nRep = 10L, nDraw = 2000L),
  Run("sim n=40 conflict m=2", MakeTarget(40, 2), nRep = 8L,  nDraw = 2500L),
  Run("sim n=40 conflict m=8", MakeTarget(40, 8), nRep = 8L,  nDraw = 2500L))
saveRDS(res, file.path(SCR, "vine-ls-probe7.rds"))
cat("\n\n==================== SUMMARY ====================\n")
print(res[, c("label", "arm", "patience", "pHit", "medEvals", "evalsSearch",
              "evalRatio", "pCensored")], row.names = FALSE, digits = 4)
