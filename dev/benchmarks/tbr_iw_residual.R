# tbr_iw_residual.R -- characterise an IW unrooted-TBR residual miss.
# For a given (nTip, treeIndex), reproduce the kernel's converged tree, then
# split the improving neighbourhood into SPR-reachable vs TBR-only, and report
# the broken bipartition of the single best improving move (vs the converged
# tree) so we can tell if it is the ROOT edge or a non-root edge, and whether
# the deficit is at the SPR (clip+graft) level or the TBR (reroot) level.
#
# Usage: Rscript dev/benchmarks/tbr_iw_residual.R [nTip] [treeIndex] [concavity]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly = TRUE)
nTip <- if (length(args) >= 1) as.integer(args[[1]]) else 16L
idx  <- if (length(args) >= 2) as.integer(args[[2]]) else 19L
conc <- if (length(args) >= 3) as.numeric(args[[3]]) else 10
nChar <- 60L; nState <- 3L; eps <- 1e-6

randomData <- function(seed) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  m <- matrix(sample(0:(nState - 1L), nTip * nChar, replace = TRUE),
              nrow = nTip, dimnames = list(tips, NULL))
  phy <- phangorn::phyDat(m, type = "USER", levels = as.character(0:(nState - 1L)))
  at <- attributes(phy)
  list(phy = phy, contrast = at$contrast,
       tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, nTip = length(phy), labels = names(phy))
}
scoreTree <- function(tree, d) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  TreeSearch:::ts_fitch_score(edge, d$contrast, d$tip_data, d$weight, d$levels, concavity = conc)
}
kernelTbr <- function(tree, d) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, concavity = conc, unrooted = TRUE)
  structure(list(edge = res$edge, Nnode = d$nTip - 1L, tip.label = d$labels), class = "phylo")
}
# Bipartitions (as sorted tip-sets of the smaller side) of an unrooted tree,
# computed manually from the edge matrix (robust to TreeTools S3 export quirks).
bips <- function(tree) {
  tree <- Preorder(tree)
  e <- tree$edge; nTip <- length(tree$tip.label); N <- nTip + tree$Nnode
  desc <- vector("list", N)
  for (i in seq_len(nTip)) desc[[i]] <- i
  for (k in rev(seq_len(nrow(e)))) {        # reverse preorder => children first
    p <- e[k, 1]; ch <- e[k, 2]
    desc[[p]] <- c(desc[[p]], desc[[ch]])
  }
  tips <- tree$tip.label; out <- character(0)
  for (k in seq_len(nrow(e))) {
    ch <- e[k, 2]
    if (ch <= nTip) next                    # leaf edge = trivial split
    a <- sort(unique(desc[[ch]]))
    if (length(a) < 2 || length(a) > nTip - 2) next
    side <- tips[a]; other <- tips[setdiff(seq_len(nTip), a)]
    s <- if (length(side) <= length(other)) side else other
    out <- c(out, paste(sort(s), collapse = ","))
  }
  sort(unique(out))
}

d <- randomData(1000L + idx)
set.seed(7000L + idx); start <- RandomTree(d$phy, root = TRUE)
set.seed(idx)
conv <- kernelTbr(start, d)
base <- scoreTree(conv, d)
cat(sprintf("=== IW residual: nTip=%d tree#%d concavity=%g ===\n", nTip, idx, conc))
cat(sprintf("converged IW = %.5f\n", base))

# Enumerate neighbours at two rootings; keep best SPR and best TBR separately.
sprScores <- c(); tbrTrees <- list(); tbrScores <- c()
for (rt in d$labels[1:2]) {
  rr <- Preorder(RootTree(conv, rt))
  sm <- SPRMoves(rr); tm <- TBRMoves(rr)
  sprScores <- c(sprScores, vapply(sm, scoreTree, double(1), d = d))
  ts <- vapply(tm, scoreTree, double(1), d = d)
  tbrScores <- c(tbrScores, ts); tbrTrees <- c(tbrTrees, tm)
}
bestSpr <- if (length(sprScores)) min(sprScores) else Inf
bestTbr <- if (length(tbrScores)) min(tbrScores) else Inf
cat(sprintf("best SPR neighbour  = %.5f  (improves: %s)\n", bestSpr, bestSpr < base - eps))
cat(sprintf("best TBR neighbour  = %.5f  (improves: %s)\n", bestTbr, bestTbr < base - eps))
cat(sprintf("=> deficit level: %s\n",
            if (bestSpr < base - eps) "SPR (basic clip+graft miss)" else
            if (bestTbr < base - eps) "TBR-reroot only" else "none (oracle artifact?)"))

# Characterise the single best improving move: which bipartition(s) does it
# add/remove vs the converged tree? Is the changed edge incident to the root?
if (bestTbr < base - eps) {
  wi <- which.min(tbrScores); bestT <- Preorder(tbrTrees[[wi]])
  cb <- bips(conv); tb <- bips(bestT)
  added <- setdiff(tb, cb); removed <- setdiff(cb, tb)
  cat(sprintf("\nbest improving tree IW = %.5f (delta %.5f)\n", tbrScores[wi], base - tbrScores[wi]))
  cat(sprintf("bipartitions changed: %d removed, %d added\n", length(removed), length(added)))
  cat("REMOVED (present in converged, gone in improver):\n")
  for (s in removed) cat("   -", s, "\n")
  cat("ADDED (new in improver):\n")
  for (s in added) cat("   +", s, "\n")
  # Re-feed the improver to the kernel: does it climb further (real basin)?
  reconv <- scoreTree(kernelTbr(bestT, d), d)
  cat(sprintf("\nkernel re-run from improver converges to IW = %.5f\n", reconv))
}
