# ---------------------------------------------------------------------------
# Vine-style embedding search for LeastSquaresTree(): feasibility probe.
#
# Question: can a Vine-style "differentiate through the NJ decoder" scheme
# improve on TreeSearch's NNI+SPR least-squares search?
#
# The decisive unknowns, in order:
#   E1  Does a gradient exist?  RSS-after-refit is (predicted) piecewise
#       CONSTANT in the input distance matrix D, because refitting discards
#       the decoder's branch lengths.  The surrogate objective -- score NJ's
#       OWN branch lengths against the target -- should be piecewise smooth.
#   E2  Cell geometry.  How does the decoded topology change as D is
#       perturbed?  Is there a perturbation scale at which NJ(D + eps) is a
#       NEAR neighbour of NJ(D), or does it jump straight to "far away"?
#   E3  Surrogate fidelity.  Does the differentiable surrogate rank topologies
#       the same way the real (refit) objective does?  If not, the gradient
#       points somewhere useless.
#   E4  Headroom.  Is there anything left to win over LeastSquaresTree()?
#
# Target matrix: Lapointe & Cucumel average consensus -- the mean path-length
# distance matrix over a set of most-parsimonious trees.  Non-additive by
# construction, which is exactly the case LeastSquaresTree() was written for.
# ---------------------------------------------------------------------------

suppressMessages({
  library(TreeSearch, lib.loc = if (dir.exists(Sys.getenv("TS_LIB", ".agent-p0"))) Sys.getenv("TS_LIB", ".agent-p0") else .libPaths())
  library(ape)
})
set.seed(1)

TS <- "C:/Users/pjjg18/GitHub/TreeSearch"
OUT <- file.path(Sys.getenv("SCRATCH", unset = "."), "vine-ls-probe.rds")

# ---- target: average-consensus distance matrix over an MPT set -------------
data(inapplicable.phyData, package = "TreeSearch")
dataset <- inapplicable.phyData[["Longrich2010"]]
nTip <- length(dataset)
cat(sprintf("Dataset: Longrich2010, %d tips\n", nTip))

mpts <- MaximizeParsimony(dataset, maxReplicates = 6L, maxSeconds = 60,
                          nThreads = 2L, verbosity = 0)
mpts <- if (inherits(mpts, "phylo")) list(mpts) else
  lapply(seq_along(mpts), function(i) mpts[[i]])
cat(sprintf("MPT set: %d trees\n", length(mpts)))

labs <- sort(mpts[[1]][["tip.label"]])
PathDist <- function(tr) {
  tr[["edge.length"]] <- rep(1, nrow(tr[["edge"]]))
  cophenetic(tr)[labs, labs]
}
Dobs <- Reduce(`+`, lapply(mpts, PathDist)) / length(mpts)
diag(Dobs) <- 0
cat(sprintf("Target D: mean %.3f, sd %.3f; additive? RSS(NJ) below\n",
            mean(Dobs[upper.tri(Dobs)]), sd(Dobs[upper.tri(Dobs)])))

# ---- scoring ---------------------------------------------------------------
pairIdx <- upper.tri(Dobs)
targetVec <- Dobs[pairIdx]

# Real objective: refit branch lengths optimally on the decoded topology.
RefitRSS <- function(tr, method = "ols") {
  fit <- try(LeastSquaresFit(tr, Dobs, method = method), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  as.numeric(attr(fit, "RSS"))
}

# Surrogate objective: score the DECODER's own branch lengths.  This is the
# only version of the objective that can have a non-zero gradient wrt D.
SurrogateRSS <- function(tr) {
  if (is.null(tr[["edge.length"]])) return(NA_real_)
  cp <- cophenetic(tr)[labs, labs]
  sum((cp[pairIdx] - targetVec)^2)
}

Decode <- function(D) {
  tr <- ape::nj(as.dist(D))
  tr[["edge.length"]][tr[["edge.length"]] < 0] <- 0  # NJ can emit negatives
  tr
}

RF <- function(a, b) as.numeric(TreeDist::RobinsonFoulds(a, b))
CID <- function(a, b) as.numeric(TreeDist::ClusteringInfoDistance(a, b,
                                                                  normalize = TRUE))

# symmetric zero-diagonal perturbation, scaled to the spread of Dobs
sdD <- sd(Dobs[pairIdx])
Perturb <- function(D, eps, delta = NULL) {
  if (is.null(delta)) {
    e <- matrix(0, nTip, nTip)
    e[pairIdx] <- rnorm(sum(pairIdx))
    delta <- e + t(e)
  }
  Dp <- D + eps * sdD * delta
  Dp[Dp < 0] <- 0
  diag(Dp) <- 0
  dimnames(Dp) <- dimnames(D)
  Dp
}

njBase <- Decode(Dobs)
rssBase <- RefitRSS(njBase)
cat(sprintf("\nBaseline A  NJ(Dobs), refit OLS      RSS = %.6f\n", rssBase))
cat(sprintf("Baseline A' NJ(Dobs), NJ's own lengths RSS = %.6f  (surrogate)\n",
            SurrogateRSS(njBase)))

tSearch <- system.time(
  lsTree <- LeastSquaresTree(Dobs, method = "ols", spr = TRUE)
)[["elapsed"]]
rssSearch <- as.numeric(attr(lsTree, "RSS"))
cat(sprintf("Baseline B  LeastSquaresTree (NNI+SPR) RSS = %.6f  [%.2f s, RF to NJ = %g]\n",
            rssSearch, tSearch, RF(njBase, lsTree)))

# ---------------------------------------------------------------------------
# E1: is the refit objective piecewise constant along a line in D-space?
# ---------------------------------------------------------------------------
cat("\n--- E1: objective along a 1-D path in D-space ---\n")
e <- matrix(0, nTip, nTip); e[pairIdx] <- rnorm(sum(pairIdx))
delta1 <- e + t(e)
tGrid <- seq(0, 0.60, length.out = 241)
e1 <- do.call(rbind, lapply(tGrid, function(t) {
  tr <- Decode(Perturb(Dobs, t, delta1))
  data.frame(t = t, refit = RefitRSS(tr), surr = SurrogateRSS(tr),
             rf = RF(njBase, tr))
}))
# count distinct plateaux in the refit series
runs <- rle(round(e1$refit, 8))
cat(sprintf("refit RSS: %d distinct values over %d grid points (%d plateaux)\n",
            length(unique(round(e1$refit, 8))), nrow(e1), length(runs$lengths)))
cat(sprintf("surrogate: %d distinct values over %d grid points\n",
            length(unique(round(e1$surr, 8))), nrow(e1)))
# largest jump between ADJACENT grid points, refit vs surrogate
cat(sprintf("max adjacent step  refit = %.4f   surrogate = %.4f\n",
            max(abs(diff(e1$refit))), max(abs(diff(e1$surr)))))
cat(sprintf("median adjacent step refit = %.6f  surrogate = %.6f\n",
            median(abs(diff(e1$refit))), median(abs(diff(e1$surr)))))

# ---------------------------------------------------------------------------
# E2: cell geometry -- topology change vs perturbation scale
# ---------------------------------------------------------------------------
cat("\n--- E2: decoded-topology change vs perturbation scale ---\n")
epsGrid <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0)
nDraw <- 120
maxSplits <- 2 * (nTip - 3)
e2 <- do.call(rbind, lapply(epsGrid, function(eps) {
  res <- do.call(rbind, lapply(seq_len(nDraw), function(i) {
    tr <- Decode(Perturb(Dobs, eps))
    data.frame(eps = eps, rf = RF(njBase, tr), cid = CID(njBase, tr),
               refit = RefitRSS(tr), surr = SurrogateRSS(tr))
  }))
  res
}))
agg <- do.call(rbind, lapply(split(e2, e2$eps), function(d) {
  data.frame(eps = d$eps[1],
             pSame = mean(d$rf == 0),
             pNear = mean(d$rf > 0 & d$rf <= 4),   # <= 2 NNI-ish moves away
             medRF = median(d$rf[d$rf > 0]),
             medRFall = median(d$rf),
             maxRF = max(d$rf),
             medCID = median(d$cid),
             bestRefit = min(d$refit, na.rm = TRUE))
}))
agg$medRF[is.na(agg$medRF)] <- 0
print(agg, row.names = FALSE, digits = 3)
cat(sprintf("(total unrooted splits = %d; RF is a count of differing splits)\n",
            maxSplits))

# ---------------------------------------------------------------------------
# E3: surrogate fidelity
# ---------------------------------------------------------------------------
cat("\n--- E3: does the surrogate rank topologies like the real objective? ---\n")
ok <- is.finite(e2$refit) & is.finite(e2$surr)
cat(sprintf("pooled Spearman(surrogate, refit) = %.3f  (n = %d)\n",
            cor(e2$surr[ok], e2$refit[ok], method = "spearman"), sum(ok)))
byEps <- do.call(rbind, lapply(split(e2[ok, ], e2$eps[ok]), function(d) {
  data.frame(eps = d$eps[1], n = nrow(d),
             spearman = if (nrow(d) > 5 && sd(d$refit) > 0)
               cor(d$surr, d$refit, method = "spearman") else NA_real_)
}))
print(byEps, row.names = FALSE, digits = 3)

# does the surrogate's argmin coincide with the refit argmin?
bestSurr <- e2[ok, ][which.min(e2$surr[ok]), ]
bestRefit <- e2[ok, ][which.min(e2$refit[ok]), ]
cat(sprintf("argmin surrogate -> refit RSS %.5f ; argmin refit -> refit RSS %.5f\n",
            bestSurr$refit, bestRefit$refit))

# ---------------------------------------------------------------------------
# E4: headroom -- random search in D-space vs NNI+SPR
# ---------------------------------------------------------------------------
cat("\n--- E4: headroom ---\n")
epsBest <- agg$eps[which.min(agg$bestRefit)]
cat(sprintf("using eps = %g (best refit RSS found in E2)\n", epsBest))
nRand <- 400
tRand <- system.time({
  randRSS <- vapply(seq_len(nRand), function(i)
    RefitRSS(Decode(Perturb(Dobs, epsBest))), numeric(1))
})[["elapsed"]]
cat(sprintf("random D-space search, %d draws: best RSS = %.6f  [%.2f s]\n",
            nRand, min(randRSS, na.rm = TRUE), tRand))
cat(sprintf("   vs NJ start        RSS = %.6f\n", rssBase))
cat(sprintf("   vs NNI+SPR search  RSS = %.6f\n", rssSearch))

# multi-start: are perturbed-NJ trees good STARTING trees for the existing search?
nStart <- 12
starts <- lapply(seq_len(nStart), function(i) Decode(Perturb(Dobs, epsBest)))
tMulti <- system.time({
  msRSS <- vapply(starts, function(s) {
    as.numeric(attr(LeastSquaresTree(Dobs, tree = s, method = "ols"), "RSS"))
  }, numeric(1))
})[["elapsed"]]
cat(sprintf("multi-start NNI+SPR from %d perturbed-NJ trees: best = %.6f, ",
            nStart, min(msRSS)))
cat(sprintf("median = %.6f  [%.2f s]\n", median(msRSS), tMulti))

saveRDS(list(Dobs = Dobs, e1 = e1, e2 = e2, agg = agg, byEps = byEps,
             rssBase = rssBase, rssSearch = rssSearch, randRSS = randRSS,
             msRSS = msRSS, epsBest = epsBest, nTip = nTip,
             nMPT = length(mpts)), OUT)
cat(sprintf("\nsaved -> %s\n", OUT))
