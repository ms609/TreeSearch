# tbr_oracle_iw.R -- differential completeness oracle that scores with the
# KERNEL's own scorer (ts_fitch_score, with concavity), so it is valid for IW
# (and EW) without any TreeLength scoring-match.  Trajectory-independent: run the
# in-kernel direct unrooted TBR to convergence, then assert no TBR/SPR neighbour
# (TreeTools enumerators, two rootings) scores strictly lower by ts_fitch_score.
#
# Usage: Rscript dev/benchmarks/tbr_oracle_iw.R [nTrees] [nTip] [concavity]
#   concavity < 0 => EW;  finite > 0 => IW
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})
args   <- commandArgs(trailingOnly = TRUE)
nTrees <- if (length(args) >= 1) as.integer(args[[1]]) else 60L
nTip   <- if (length(args) >= 2) as.integer(args[[2]]) else 12L
conc   <- if (length(args) >= 3) as.numeric(args[[3]]) else -1
unroot <- if (length(args) >= 4) as.logical(as.integer(args[[4]])) else TRUE
nChar  <- 60L; nState <- 3L
eps    <- 1e-6

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

# Kernel-matched score of any tree (EW or IW, by concavity).
scoreTree <- function(tree, d) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  TreeSearch:::ts_fitch_score(edge, d$contrast, d$tip_data, d$weight, d$levels,
                              concavity = conc)
}

kernelTbr <- function(tree, d) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L,
    concavity = conc, unrooted = unroot)
  structure(list(edge = res$edge, Nnode = d$nTip - 1L, tip.label = d$labels),
            class = "phylo")
}

# Best improving neighbour by ts_fitch_score; all_tbr/all_spr at two rootings
# cover every break edge.  Returns improving score, or NULL if clean.
bestImproving <- function(tree, d) {
  base <- scoreTree(tree, d)
  cand <- list()
  for (rt in d$labels[1:2]) {
    rr <- Preorder(RootTree(tree, rt))
    cand <- c(cand, TBRMoves(rr), SPRMoves(rr))
  }
  if (!length(cand)) return(NULL)
  ls <- vapply(cand, scoreTree, double(1), d = d)
  if (min(ls) < base - eps) list(len = min(ls), base = base) else NULL
}

cat(sprintf("=== IW/EW differential oracle (ts_fitch_score): %d trees, %d tips, concavity=%s (%s) ===\n",
            nTrees, nTip, conc, if (conc < 0) "EW" else "IW"))
fails <- 0L; firstFail <- NULL
for (i in seq_len(nTrees)) {
  d <- randomData(1000L + i)
  set.seed(7000L + i); start <- RandomTree(d$phy, root = TRUE)
  set.seed(i)
  conv <- kernelTbr(start, d)
  imp <- bestImproving(conv, d)
  if (!is.null(imp)) {
    fails <- fails + 1L
    if (is.null(firstFail)) firstFail <- list(i = i, conv = scoreTree(conv, d), imp = imp)
  }
}
cat(sprintf("\nRESULT: %d / %d converged trees had an improving neighbour (FAILURES)\n", fails, nTrees))
if (fails == 0L) cat("=> DIRECT unrooted TBR is COMPLETE on this sample (kernel-scored).\n") else {
  ff <- firstFail
  cat(sprintf("=> INCOMPLETE. First: tree #%d, converged %.4f, neighbour %.4f (miss %.4f)\n",
              ff$i, ff$conv, ff$imp$len, ff$conv - ff$imp$len))
}
