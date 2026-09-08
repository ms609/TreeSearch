# tbr_collapsed_test.R -- single-variable test (advisor, 2026-06-18).
#
# Confirms collapsed-edge pruning is the residual NEIGHBOURHOOD cause: it drops
# strict-improving TBR moves because its "provably cannot improve" proof is the
# SPR one and is unsound for TBR fragment-rerooting.
#
# TBRParams::unrooted=TRUE currently DISABLES collapsed pruning (only).  The
# all-tips reroot emulation supplies the rooting coverage, so this isolates
# collapsed as the one variable.  SUCCESS CRITERION: the resulting optimum's
# canonical-TBR neighbourhood has 0 improving neighbours (enumerator-clean).
source("dev/benchmarks/tbr_shared_start_lib.R")
d <- prepareDataset("Zanol2014")
norm    <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))

# TsTbr with collapsed pruning OFF (unrooted=TRUE).
TsTbrU <- function(startTree, seed, acceptEqual = FALSE, maxHits = 1L) {
  edge <- PhyloToKernelEdge(startTree, d); set.seed(seed)
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = maxHits, acceptEqual = acceptEqual, maxChanges = 0L,
    unrooted = TRUE)
  tr <- structure(list(edge = res$edge, Nnode = d$nTip - 1L,
                        tip.label = names(d$phy)), class = "phylo")
  list(tree = norm(tr), len = TreeLength(tr, d$phy))
}

# All-tips reroot-invariant TBR, collapsed OFF.
RootInvariantU <- function(startTree, seed, rerootTips = names(d$phy)) {
  cur <- TsTbrU(startTree, seed); best <- cur$tree; bestLen <- cur$len
  repeat {
    improved <- FALSE
    for (tp in rerootTips) {
      rr <- norm(ape::root(best, outgroup = tp, resolve.root = TRUE))
      r  <- TsTbrU(rr, seed)
      if (r$len < bestLen) { bestLen <- r$len; best <- r$tree; improved <- TRUE }
    }
    if (!improved) break
  }
  list(len = bestLen, tree = best)
}

probe0 <- function(tree, label, baseLen) {
  nb <- TBRMoves(norm(tree)); ls <- vapply(nb, TreeLength, double(1), d$phy)
  best <- min(ls); nBetter <- sum(ls < baseLen - 0.5)
  cat(sprintf("  [%s] enumerator: %d neighbours, best=%.0f, %d improving %s\n",
              label, length(nb), best, nBetter,
              if (nBetter == 0) "=> CLEAN (0 improving)" else "=> still incomplete"))
  invisible(list(best = best, nBetter = nBetter))
}

cat("=== Collapsed-pruning isolation test (Zanol2014) ===\n\n")

# (A) Direct: feed the OLD collapsed-ON optimum (1284) into collapsed-OFF descent.
optC <- norm(ape::read.tree("dev/benchmarks/tbr_results/ts_reroot_invariant_opt.tre"))
cat(sprintf("(A) collapsed-ON optimum = %.0f\n", TreeLength(optC, d$phy)))
a1 <- TsTbrU(optC, seed = 1)
cat(sprintf("    -> collapsed-OFF strict descent (same rooting): %.0f -> %.0f\n",
            TreeLength(optC, d$phy), a1$len))
aRI <- RootInvariantU(optC, seed = 1)
cat(sprintf("    -> collapsed-OFF all-tips reroot-invariant:      %.0f -> %.0f\n",
            TreeLength(optC, d$phy), aRI$len))
probe0(aRI$tree, "A all-tips opt", aRI$len)

# (B) From scratch: random seed1 start, collapsed-OFF all-tips reroot-invariant.
st <- norm({ set.seed(1001); RandomTree(d$phy, root = TRUE) })
cat(sprintf("\n(B) random start = %.0f\n", TreeLength(st, d$phy)))
bRI <- RootInvariantU(st, seed = 1)
cat(sprintf("    -> collapsed-OFF all-tips reroot-invariant: %.0f\n", bRI$len))
probe0(bRI$tree, "B all-tips opt", bRI$len)
ape::write.tree(bRI$tree, "dev/benchmarks/tbr_results/ts_collapsedoff_opt.tre")

cat("\nReading: optimum dropping below 1284 AND probe 0-improving => collapsed pruning\n",
    "was the residual neighbourhood cause (SPR-sound, TBR-unsound). Remaining gap to\n",
    "TNT's 1264 (if any) is then basin/path, not move-set.\n", sep = "")
