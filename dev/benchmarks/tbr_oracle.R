# tbr_oracle.R -- small-tree DIFFERENTIAL ORACLE for kernel TBR completeness.
#
# "Everything correct" (project lead, 2026-06-18) made checkable: run the
# IN-KERNEL tbr_search (the path the default actually uses -- single rooting, no
# R reroot scaffold) to convergence, then assert the package's own UNOPTIMISED
# enumerator (TBRMoves/SPRMoves -> all_tbr/all_spr in rearrange.cpp) finds NO
# strictly-improving neighbour.  Any failure prints a concrete small tree + the
# exact missing move -- debuggable in seconds, and a permanent regression test.
#
# Scope: tests the cpp-search kernel (ts_tbr_diagnostics -> tbr_search), which is
# what MaximizeParsimony / driven / ratchet / sector / fuse all route through.
# The old CustomSearch/TreeSearch path (RootedTBRSwap + Morphy) is separate.
#
# Usage: Rscript dev/benchmarks/tbr_oracle.R [nTrees=200] [nTip=12] [unrooted=0]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})

args     <- commandArgs(trailingOnly = TRUE)
nTrees   <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
nTip     <- if (length(args) >= 2) as.integer(args[[2]]) else 12L
unrooted <- if (length(args) >= 3) as.logical(as.integer(args[[3]])) else FALSE
emul     <- if (length(args) >= 4) as.logical(as.integer(args[[4]])) else FALSE
nChar    <- 60L
nState   <- 3L          # states 0,1,2 (+ '?')

# Build a random USER-type phyDat on `nTip` tips, `nChar` characters.
randomData <- function(seed) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  m <- matrix(sample(0:(nState - 1L), nTip * nChar, replace = TRUE),
              nrow = nTip, dimnames = list(tips, NULL))
  phy <- phangorn::phyDat(m, type = "USER", levels = as.character(0:(nState - 1L)))
  at <- attributes(phy)
  list(phy = phy, contrast = at$contrast,
       tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, nTip = length(phy),
       labels = names(phy))
}

# Run the in-kernel TBR (single rooting -- as the default runs) to convergence.
# `unrooted=TRUE` requires the move-fix prototype build; on the clean post-fix
# cpp-search build (no `unrooted` arg) we call the default signature.
kernelTbr <- function(tree, d, unrooted) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  res <- if (isTRUE(unrooted))
    TreeSearch:::ts_tbr_diagnostics(
      edge, d$contrast, d$tip_data, d$weight, d$levels,
      maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, unrooted = TRUE)
  else
    TreeSearch:::ts_tbr_diagnostics(
      edge, d$contrast, d$tip_data, d$weight, d$levels,
      maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L)
  structure(list(edge = res$edge, Nnode = d$nTip - 1L, tip.label = d$labels),
            class = "phylo")
}

# All-tips reroot emulation: converge, then reroot at each tip and re-converge,
# adopting improvements, until a full tip sweep yields nothing.  Approximates a
# true in-kernel unrooted/root-edge mechanism using the single-rooting kernel.
# Isolates the ROOT-EDGE residual from any move-generation residual.
kernelTbrEmul <- function(tree, d, unrooted) {
  best <- kernelTbr(tree, d, unrooted); bestLen <- TreeLength(best, d$phy)
  repeat {
    improved <- FALSE
    for (tp in d$labels) {
      rr <- Preorder(RootTree(best, tp))
      r  <- kernelTbr(rr, d, unrooted); rl <- TreeLength(r, d$phy)
      if (rl < bestLen - 0.5) { bestLen <- rl; best <- r; improved <- TRUE }
    }
    if (!improved) break
  }
  best
}

# Full unrooted-TBR cleanliness check: all_tbr at TWO distinct rootings (tip1 &
# tip2) covers every break edge (each rooting only omits its own root-edge =
# that tip's pendant); plus all_spr for good measure.  Returns the best
# improving neighbour length and tree, or NULL if clean.
bestImproving <- function(tree, d) {
  base <- TreeLength(tree, d$phy)
  cand <- list()
  for (rt in d$labels[1:2]) {
    rr <- Preorder(RootTree(tree, rt))
    cand <- c(cand, TBRMoves(rr), SPRMoves(rr))
  }
  if (!length(cand)) return(NULL)
  ls <- vapply(cand, TreeLength, double(1), d$phy)
  if (min(ls) < base - 0.5)
    list(len = min(ls), tree = cand[[which.min(ls)]], base = base)
  else NULL
}

cat(sprintf("=== TBR completeness oracle: %d trees, %d tips, unrooted=%s ===\n",
            nTrees, nTip, unrooted))
fails <- 0L; firstFail <- NULL
for (i in seq_len(nTrees)) {
  d <- randomData(1000L + i)
  set.seed(7000L + i)
  start <- RandomTree(d$phy, root = TRUE)
  set.seed(i)                       # kernel RNG (clip order)
  conv <- if (emul) kernelTbrEmul(start, d, unrooted) else kernelTbr(start, d, unrooted)
  convLen <- TreeLength(conv, d$phy)
  imp <- bestImproving(conv, d)
  if (!is.null(imp)) {
    fails <- fails + 1L
    if (is.null(firstFail))
      firstFail <- list(i = i, conv = conv, convLen = convLen, imp = imp,
                        start = start, d = d)
  }
}
cat(sprintf("\nRESULT: %d / %d converged trees had an improving canonical neighbour (FAILURES)\n",
            fails, nTrees))
if (fails == 0L) {
  cat("=> KERNEL TBR IS COMPLETE on this sample (0 missing moves).\n")
} else {
  ff <- firstFail
  cat(sprintf("=> INCOMPLETE. First failure: tree #%d, kernel converged %.0f, canonical finds %.0f (miss %.0f)\n",
              ff$i, ff$convLen, ff$imp$len, ff$convLen - ff$imp$len))
  cat("   converged Newick:", ape::write.tree(Preorder(RootTree(ff$conv, ff$d$labels[1]))), "\n")
  cat("   improving  Newick:", ape::write.tree(Preorder(RootTree(ff$imp$tree, ff$d$labels[1]))), "\n")
  saveRDS(ff, "dev/benchmarks/tbr_results/oracle_first_fail.rds")
  cat("   (saved first failure -> dev/benchmarks/tbr_results/oracle_first_fail.rds)\n")
}
