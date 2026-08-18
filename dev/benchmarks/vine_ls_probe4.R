# ---------------------------------------------------------------------------
# Probe 4 -- final, on a REPRODUCIBLE and TUNABLY HARD target family.
#
# Probes 2/3 were confounded: the average-consensus target was built from a
# stochastic MaximizeParsimony() run, so the MPT set (and hence the problem's
# difficulty) changed between probes -- Wortley2006 gave 19 MPTs in probe 2 and
# 2 in probe 3.  Worse, MPTs of a single dataset are all similar, so their
# average path-length matrix is nearly additive and NJ almost solves it.
#
# Here the target is built from a controlled set of conflicting trees:
#   backbone T0, plus (k-1) copies each perturbed by `m` random SPR moves.
#   Dobs = mean path-length matrix.   Conflict level is set by `m`.
# Fully reproducible from a seed; difficulty is a dial, not an accident.
# One cached real MPT-set target is retained for external validity.
#
# Per target, all under NNLS (the default and the Lapointe & Cucumel convention):
#   [1] NJ(Dobs)                                        -- the free baseline
#   [2] LeastSquaresTree(): NNI+SPR                     -- the incumbent
#   [3] random D-space search, TIME-MATCHED to [2]      -- undirected Vine move
#   [4] ORACLE greedy D-space hill-climb on the TRUE objective, time-matched
#         == the best any gradient could do, since a perfect gradient only
#            tells you which way the true objective improves.
#   [5] surrogate ranking among DISTINCT immediate (RF==2) neighbours
#         == whether the only differentiable objective can pick that direction
#   [6] OLS-vs-NNLS cross-scoring of each winner
# ---------------------------------------------------------------------------

suppressMessages({
  library(TreeSearch, lib.loc = if (dir.exists(Sys.getenv("TS_LIB", ".agent-p0"))) Sys.getenv("TS_LIB", ".agent-p0") else .libPaths())
  library(ape)
})
SCR <- Sys.getenv("SCRATCH", unset = ".")
OUT <- file.path(SCR, "vine-ls-probe4.rds")
METHOD <- "nnls"

# ---- reproducible target construction --------------------------------------
MakeTarget <- function(nTip, m, k = 12L, seed = 7L) {
  set.seed(seed)
  t0 <- TreeTools::RandomTree(nTip, root = FALSE)
  t0[["tip.label"]] <- sprintf("t%02d", seq_len(nTip))
  trees <- c(list(t0), lapply(seq_len(k - 1L), function(i)
    phangorn::rSPR(t0, moves = m)))
  labs <- sort(t0[["tip.label"]])
  P <- lapply(trees, function(tr) {
    tr[["edge.length"]] <- rep(1, nrow(tr[["edge"]])); cophenetic(tr)[labs, labs]
  })
  D <- Reduce(`+`, P) / length(P); diag(D) <- 0
  list(D = D, labs = labs, trees = trees)
}

CachedMPTTarget <- function(dsName, seed = 11L) {
  f <- file.path(SCR, paste0("mpt-", dsName, ".rds"))
  if (file.exists(f)) { mpts <- readRDS(f) } else {
    data(inapplicable.phyData, package = "TreeSearch")
    set.seed(seed)
    r <- MaximizeParsimony(inapplicable.phyData[[dsName]], maxReplicates = 8L,
                           maxSeconds = 120, nThreads = 2L, verbosity = 0)
    mpts <- if (inherits(r, "phylo")) list(r) else
      lapply(seq_along(r), function(i) r[[i]])
    saveRDS(mpts, f)
  }
  labs <- sort(mpts[[1]][["tip.label"]])
  P <- lapply(mpts, function(tr) {
    tr[["edge.length"]] <- rep(1, nrow(tr[["edge"]])); cophenetic(tr)[labs, labs]
  })
  D <- Reduce(`+`, P) / length(P); diag(D) <- 0
  list(D = D, labs = labs, trees = mpts)
}

# ---- one full evaluation ----------------------------------------------------
Evaluate <- function(label, tgt, eps = 0.15, nNearWanted = 250L) {
  D <- tgt$D; labs <- tgt$labs; nTip <- length(labs)
  pairIdx <- upper.tri(D); targetVec <- D[pairIdx]; sdD <- sd(targetVec)
  cat(sprintf("\n=========== %s | %d tips, %d source trees ===========\n",
              label, nTip, length(tgt$trees)))

  Decode <- function(M) {
    tr <- ape::nj(as.dist(M)); tr[["edge.length"]][tr[["edge.length"]] < 0] <- 0; tr
  }
  Perturb <- function(M, e) {
    z <- matrix(0, nTip, nTip); z[pairIdx] <- rnorm(sum(pairIdx))
    Mp <- M + e * sdD * (z + t(z)); Mp[Mp < 0] <- 0; diag(Mp) <- 0
    dimnames(Mp) <- dimnames(D); Mp
  }
  Refit <- function(tr, meth = METHOD) {
    f <- try(LeastSquaresFit(tr, D, method = meth), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(attr(f, "RSS"))
  }
  Surr <- function(tr) sum((cophenetic(tr)[labs, labs][pairIdx] - targetVec)^2)
  RF <- function(a, b) as.numeric(TreeDist::RobinsonFoulds(a, b))

  nj <- Decode(D); rssNJ <- Refit(nj)
  tS <- system.time(best <- LeastSquaresTree(D, method = METHOD))[["elapsed"]]
  rssS <- as.numeric(attr(best, "RSS"))
  cat(sprintf("[1] NJ(Dobs)              RSS = %10.4f\n", rssNJ))
  cat(sprintf("[2] NNI+SPR search        RSS = %10.4f   [%.1f s, RF to NJ = %g]\n",
              rssS, tS, RF(nj, best)))
  budget <- max(tS, 2)

  # [3] time-matched undirected D-space search
  t0 <- proc.time()[["elapsed"]]; bR <- rssNJ; nR <- 0L
  repeat {
    v <- Refit(Decode(Perturb(D, eps))); nR <- nR + 1L
    if (is.finite(v) && v < bR) bR <- v
    if (proc.time()[["elapsed"]] - t0 > budget) break
  }
  cat(sprintf("[3] random D-space        RSS = %10.4f   [%d draws in %.1f s]\n",
              bR, nR, budget))

  # [4] ORACLE greedy hill-climb in D-space on the TRUE objective
  Dc <- D; bO <- rssNJ; t0 <- proc.time()[["elapsed"]]
  nP <- 0L; nA <- 0L; stall <- 0L; e <- eps
  repeat {
    Dp <- Perturb(Dc, e); v <- Refit(Decode(Dp)); nP <- nP + 1L
    if (is.finite(v) && v < bO - 1e-12) { bO <- v; Dc <- Dp; nA <- nA + 1L; stall <- 0L }
    else { stall <- stall + 1L
           if (stall > 50L) { e <- e * 1.5; stall <- 0L }
           if (e > 1.2) { e <- eps; Dc <- D } }
    if (proc.time()[["elapsed"]] - t0 > budget) break
  }
  verdict <- if (bO < rssS - 1e-9) "ORACLE WINS" else
    if (abs(bO - rssS) < 1e-9) "oracle ties" else "ORACLE LOSES"
  cat(sprintf("[4] ORACLE D-space climb  RSS = %10.4f   [%d proposals, %d accepted]  -> %s\n",
              bO, nP, nA, verdict))

  # [5] surrogate ranking among DISTINCT immediate neighbours
  keys <- character(0); rec <- list(); tries <- 0L
  while (length(rec) < nNearWanted && tries < nNearWanted * 40L) {
    tries <- tries + 1L
    tr <- Decode(Perturb(D, eps))
    if (RF(nj, tr) != 2) next
    key <- paste(sort(as.character(TreeTools::as.Splits(tr))), collapse = "|")
    if (key %in% keys) next
    keys <- c(keys, key)
    rec[[length(rec) + 1L]] <- data.frame(refit = Refit(tr), surr = Surr(tr))
  }
  nb <- do.call(rbind, rec); nb <- nb[is.finite(nb$refit), ]
  srho <- NA_real_; pickTop <- NA; pickV <- NA_real_
  if (!is.null(nb) && nrow(nb) > 10) {
    srho <- cor(nb$surr, nb$refit, method = "spearman")
    pickV <- nb$refit[which.min(nb$surr)]
    top <- sort(nb$refit)[max(1, ceiling(0.1 * nrow(nb)))]
    pickTop <- pickV <= top
    cat(sprintf("[5] surrogate on %3d distinct RF==2 neighbours: Spearman %+.3f | picks %.4f (best %.4f) | top decile: %s | %.0f%% beat NJ\n",
                nrow(nb), srho, pickV, min(nb$refit),
                if (pickTop) "YES" else "NO", 100 * mean(nb$refit < rssNJ)))
  } else cat("[5] too few distinct RF==2 neighbours\n")

  # [6] OLS vs NNLS cross-scoring
  bestOLS <- LeastSquaresTree(D, method = "ols")
  x <- data.frame(tree = c("NJ", "NNLS winner", "OLS winner"),
                  ols = c(Refit(nj, "ols"), Refit(best, "ols"), Refit(bestOLS, "ols")),
                  nnls = c(Refit(nj, "nnls"), Refit(best, "nnls"), Refit(bestOLS, "nnls")))
  cat("[6] cross-scored:\n"); print(x, row.names = FALSE, digits = 6)

  data.frame(label = label, nTip = nTip, nSrc = length(tgt$trees),
             NJ = rssNJ, NNISPR = rssS, randD = bR, oracleD = bO,
             tSearch = tS, oracleVerdict = verdict, surrRho = srho,
             surrPick = pickV, surrTopDecile = pickTop,
             olsWinner_nnls = x$nnls[3], nnlsWinner_nnls = x$nnls[2])
}

set.seed(101)
rows <- list()
rows[[1]] <- Evaluate("sim n=20 conflict m=2", MakeTarget(20, m = 2))
rows[[2]] <- Evaluate("sim n=20 conflict m=8", MakeTarget(20, m = 8))
rows[[3]] <- Evaluate("sim n=40 conflict m=2", MakeTarget(40, m = 2))
rows[[4]] <- Evaluate("sim n=40 conflict m=8", MakeTarget(40, m = 8))
rows[[5]] <- Evaluate("real Longrich2010 MPTs", CachedMPTTarget("Longrich2010"))

S <- do.call(rbind, rows)
saveRDS(S, OUT)
cat("\n\n======================== SUMMARY ========================\n")
print(S[, c("label", "nTip", "NJ", "NNISPR", "randD", "oracleD", "oracleVerdict")],
      row.names = FALSE, digits = 6)
cat("\n-- surrogate (the only differentiable objective) --\n")
print(S[, c("label", "surrRho", "surrPick", "surrTopDecile")],
      row.names = FALSE, digits = 4)
cat("\n-- OLS-mode winner, rescored under NNLS, vs the NNLS winner --\n")
print(S[, c("label", "nnlsWinner_nnls", "olsWinner_nnls")],
      row.names = FALSE, digits = 6)
cat(sprintf("\nsaved -> %s\n", OUT))
