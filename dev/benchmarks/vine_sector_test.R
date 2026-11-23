# ---------------------------------------------------------------------------
# Sector test: is a perturbed-NJ decode a better START-TREE GENERATOR than
# random-addition-sequence (RAS) Wagner for small PARSIMONY subproblems?
#
# Everything measured in probes 1-7 was on a LEAST-SQUARES objective, where the
# decoder and the objective are natively aligned (NJ is a distance method; the
# objective is distance fit).  Sectorial search scores PARSIMONY.  None of the
# earlier numbers carry over -- this measures the parsimony case directly.
#
# Insertion point being modelled: build_ras_sector() (src/ts_sector.cpp:884),
# called in the `for (s = 0; s < ras_starts; ++s)` loop at :1186.  Sector
# starts 1..N are currently RAS Wagner trees.  pNJ competes with RAS there --
# NOT with TBR, which runs afterwards either way.
#
# Two questions, because they can disagree:
#   Q1 QUALITY   -- better sector solutions after a fixed TBR budget?
#   Q2 DIVERSITY -- RAS's randomness is the POINT; multiple starts want
#                   different basins.  pNJ decodes concentrate near the
#                   distance-optimal region, so a generator that is
#                   individually better but collectively narrower could be a
#                   net loss across restarts.
#
# Simplification, stated honestly: a real RSS sector carries the rest of the
# tree as HTU composite terminals.  Here a sector is the reduced dataset on a
# clade's taxa alone.  That is the "solve each sub-problem independently on a
# reduced dataset" core, minus the HTU anchor.
# ---------------------------------------------------------------------------

suppressMessages({
  library(TreeSearch, lib.loc = if (dir.exists(Sys.getenv("TS_LIB", ".agent-p0"))) Sys.getenv("TS_LIB", ".agent-p0") else .libPaths())
  library(ape)
})
SCR <- Sys.getenv("SCRATCH", unset = ".")
set.seed(4242)
data(inapplicable.phyData, package = "TreeSearch")

K        <- 5L      # restarts per generator (TNT uses 3; ras_starts default 1)
MAXITER  <- 100L    # fixed TBR budget, identical for every generator
MAXHITS  <- 5L
CONC     <- Inf     # equal weights

Root1 <- function(tr) TreeTools::RootTree(tr, tr[["tip.label"]][[1]])
Score <- function(tr, d) TreeLength(Root1(tr), d, concavity = CONC)

PNJStart <- function(sub, hamm, eps) {
  n <- attr(hamm, "Size"); M <- as.matrix(hamm)
  sdD <- sd(M[upper.tri(M)])
  z <- matrix(0, n, n); z[upper.tri(z)] <- rnorm(n * (n - 1) / 2)
  Mp <- M + eps * sdD * (z + t(z)); Mp[Mp < 0] <- 0; diag(Mp) <- 0
  dimnames(Mp) <- dimnames(M)
  Root1(ape::nj(as.dist(Mp)))
}
RASStart <- function(sub) {
  AdditionTree(sub, concavity = CONC, sequence = sample(names(sub)))
}
Polish <- function(start, prepped) {
  suppressWarnings(suppressMessages(
    TreeSearch(Root1(start), prepped, EdgeSwapper = TBRSwap,
               maxIter = MAXITER, maxHits = MAXHITS, verbosity = 0L)))
}
MeanCID <- function(trees) {
  if (length(trees) < 2) return(NA_real_)
  cl <- structure(trees, class = "multiPhylo")
  mean(as.numeric(TreeDist::ClusteringInfoDistance(cl, normalize = TRUE)))
}

# ---- carve sectors out of real datasets ------------------------------------
GetSectors <- function(dsName, nSector = 4L, lo = 8L, hi = 40L) {
  d <- inapplicable.phyData[[dsName]]
  ref <- MaximizeParsimony(d, maxReplicates = 4L, maxSeconds = 45,
                           nThreads = 2L, verbosity = 0)
  ref <- if (inherits(ref, "phylo")) ref else ref[[1]]
  ref <- Root1(ref)
  nT <- length(ref[["tip.label"]])
  nodes <- (nT + 1L):max(ref[["edge"]])
  sizes <- vapply(nodes, function(nd)
    length(ape::extract.clade(ref, nd)[["tip.label"]]), integer(1))
  ok <- nodes[sizes >= lo & sizes <= hi]
  if (!length(ok)) return(list())
  ok <- sample(ok, min(nSector, length(ok)))
  lapply(ok, function(nd) {
    tips <- ape::extract.clade(ref, nd)[["tip.label"]]
    list(ds = dsName, node = nd, tips = tips, sub = d[tips])
  })
}

Evaluate <- function(sec) {
  sub <- sec$sub; n <- length(sub)
  prepped <- PrepareData(sub)
  hamm <- TreeTools::Hamming(sub)
  # reference: what a thorough search finds for this sector
  refT <- MaximizeParsimony(sub, maxReplicates = 6L, maxSeconds = 30,
                            nThreads = 2L, verbosity = 0)
  refT <- if (inherits(refT, "phylo")) refT else refT[[1]]
  bestKnown <- Score(refT, sub)
  njPlain <- Score(Root1(ape::nj(hamm)), sub)

  arms <- list(
    RAS        = function() RASStart(sub),
    `pNJ .15`  = function() PNJStart(sub, hamm, 0.15),
    `pNJ rand` = function() PNJStart(sub, hamm, runif(1, 0.05, 0.30)))

  out <- lapply(names(arms), function(a) {
    starts <- lapply(seq_len(K), function(i) arms[[a]]())
    sSc <- vapply(starts, Score, numeric(1), d = sub)
    posts <- lapply(starts, Polish, prepped = prepped)
    posts <- lapply(posts, function(p) if (inherits(p, "phylo")) p else p[[1]])
    pSc <- vapply(posts, Score, numeric(1), d = sub)
    data.frame(ds = sec$ds, n = n, arm = a, bestKnown = bestKnown,
               njPlain = njPlain,
               startMed = median(sSc), startBest = min(sSc),
               postMed = median(pSc), postBest = min(pSc),
               hitBest = min(pSc) <= bestKnown + 1e-9,
               divStart = MeanCID(starts), divPost = MeanCID(posts))
  })
  do.call(rbind, out)
}

secs <- unlist(lapply(c("Sansom2010", "Wortley2006", "OLeary1999", "Griswold1999"),
                      function(x) GetSectors(x)), recursive = FALSE)
cat(sprintf("carved %d sectors (sizes: %s)\n", length(secs),
            paste(vapply(secs, function(s) length(s$sub), integer(1)), collapse = ", ")))

parts <- lapply(seq_along(secs), function(i) {
  cat(sprintf("  sector %d/%d (%s, n=%d) ... ", i, length(secs), secs[[i]]$ds,
              length(secs[[i]]$sub)))
  r <- tryCatch(Evaluate(secs[[i]]),
                error = function(e) { cat("FAILED: ", conditionMessage(e), "\n");
                                      NULL })
  if (!is.null(r)) cat("ok\n")
  r
})
nFail <- sum(vapply(parts, is.null, logical(1)))
res <- do.call(rbind, Filter(Negate(is.null), parts))
cat(sprintf("\n%d/%d sectors evaluated (%d failed)\n",
            length(secs) - nFail, length(secs), nFail))
saveRDS(res, file.path(SCR, "vine-sector-test.rds"))

cat("\n\n==================== PER-SECTOR ====================\n")
print(res, row.names = FALSE, digits = 5)

cat("\n==================== SUMMARY BY ARM ====================\n")
agg <- do.call(rbind, lapply(split(res, res$arm), function(d) data.frame(
  arm = d$arm[1], nSectors = nrow(d),
  medStartExcess = median(d$startMed - d$bestKnown),   # generator quality, raw
  medPostExcess  = median(d$postMed  - d$bestKnown),   # after fixed TBR budget
  medBestOfK     = median(d$postBest - d$bestKnown),   # what a sector keeps
  pHitBest       = mean(d$hitBest),                    # reached sector optimum
  medDivStart    = median(d$divStart, na.rm = TRUE),   # Q2: spread of starts
  medDivPost     = median(d$divPost,  na.rm = TRUE))))
print(agg, row.names = FALSE, digits = 4)
cat("\nExcess = score - best known for that sector (0 = optimal). Lower is better.\n")
cat("div* = mean pairwise normalized clustering-information distance across the K restarts.\n")
cat(sprintf("plain unperturbed NJ, median excess: %.3f\n",
            median(res$njPlain - res$bestKnown)))

# head-to-head, paired by sector
w <- reshape(res[, c("ds", "n", "arm", "postBest")], idvar = c("ds", "n"),
             timevar = "arm", direction = "wide")
names(w) <- sub("postBest\\.", "", names(w))
cat("\n==================== PAIRED best-of-K, by sector ====================\n")
print(w, row.names = FALSE, digits = 5)
for (a in c("pNJ .15", "pNJ rand")) {
  d <- w[[a]] - w[["RAS"]]
  cat(sprintf("%-9s vs RAS : better %d, tied %d, worse %d  (median diff %+.2f steps)\n",
              a, sum(d < 0), sum(d == 0), sum(d > 0), median(d)))
}
