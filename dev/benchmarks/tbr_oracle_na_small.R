# Fast small-tree NA completeness oracle.  Random synthetic data WITH
# inapplicables (so the kernel takes the has_na convergence path), small trees so
# the brute-force bestImproving() is cheap -> many starts in seconds.  For each
# start: run the in-kernel TBR (ts_tbr_diagnostics) to convergence, then assert
# TreeTools' own TBR+SPR enumerators (two rootings) find no strictly-improving
# neighbour under the kernel's own NA scorer.  Complements tbr_oracle_na.R (real
# 74/88-tip data, slow) with a fast, high-N regression signal for the root-edge
# completeness fix.
#
# Usage: Rscript dev/benchmarks/tbr_oracle_na_small.R [nStart=100] [nTip=12] [nChar=40]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})
a <- commandArgs(trailingOnly = TRUE)
nStart <- if (length(a) >= 1) as.integer(a[[1]]) else 100L
nTip   <- if (length(a) >= 2) as.integer(a[[2]]) else 12L
nChar  <- if (length(a) >= 3) as.integer(a[[3]]) else 40L
eps <- 1e-6

# Random data over states {-,0,1,2}; "-" (~25%) is inapplicable -> has_na path.
randomNaData <- function(seed) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  toks <- c("-", "0", "1", "2")
  m <- matrix(sample(toks, nTip * nChar, replace = TRUE, prob = c(.25, .25, .25, .25)),
              nrow = nTip, dimnames = list(tips, NULL))
  # Guard: drop all-inapplicable / constant columns MatrixToPhyDat may choke on.
  keep <- apply(m, 2, function(col) length(unique(col[col != "-"])) >= 1)
  m <- m[, keep, drop = FALSE]
  MatrixToPhyDat(m)
}

scoreTree <- function(tree, d) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  TreeSearch:::ts_fitch_score(edge, d$contrast, d$tip_data, d$weight, d$levels,
                              concavity = -1)
}
kernelTbr <- function(tree, d) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  res <- TreeSearch:::ts_tbr_diagnostics(edge, d$contrast, d$tip_data, d$weight,
           d$levels, maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L,
           concavity = -1, unrooted = TRUE)
  structure(list(edge = res$edge, Nnode = length(d$labels) - 1L,
                 tip.label = d$labels), class = "phylo")
}
bestImproving <- function(tree, d) {
  base <- scoreTree(tree, d); cand <- list()
  for (rt in d$labels[1:2]) {
    rr <- Preorder(RootTree(tree, rt)); cand <- c(cand, TBRMoves(rr), SPRMoves(rr))
  }
  ls <- vapply(cand, scoreTree, double(1), d = d)
  if (min(ls) < base - eps) list(len = min(ls), base = base) else NULL
}

cat(sprintf("=== small NA oracle: %d starts, %d tips, %d chars ===\n",
            nStart, nTip, nChar))
fails <- 0L; ff <- NULL
for (i in seq_len(nStart)) {
  dataset <- randomNaData(20000L + i)
  at <- attributes(dataset)
  d <- list(contrast = at$contrast,
            tip_data = matrix(unlist(dataset, use.names = FALSE),
                              nrow = length(dataset), byrow = TRUE),
            weight = at$weight, levels = at$levels, labels = names(dataset))
  set.seed(7000L + i); start <- RandomTree(dataset, root = TRUE)
  set.seed(i); conv <- kernelTbr(start, d)
  imp <- bestImproving(conv, d)
  if (!is.null(imp)) {
    fails <- fails + 1L
    if (is.null(ff)) ff <- list(i = i, conv = scoreTree(conv, d), imp = imp)
  }
}
cat(sprintf("RESULT: %d / %d converged trees had an improving neighbour\n",
            fails, nStart))
if (fails == 0L) cat("=> COMPLETE on this sample.\n") else
  cat(sprintf("=> INCOMPLETE. First: start #%d, converged %.4f, neighbour %.4f (miss %.4f)\n",
              ff$i, ff$conv, ff$imp$len, ff$conv - ff$imp$len))
