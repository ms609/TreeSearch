# tbr_neighbourhood_probe.R -- fork-settler (advisor, 2026-06-18).
#
# The reroot-invariant cross-feed showed TNT nomulpars improves the all-tips
# reroot-invariant TS optimum (1284 -> 1270).  Two possible causes:
#   (i)  TS's KERNEL TBR is incomplete: its optimisations (L812 smaller-subtree
#        skip, collapsed-edge pruning, indirect-scoring cutoff) prune real moves,
#        so it declares convergence with an improving canonical-TBR move present.
#   (ii) TNT's "TBR" exceeds textbook single-tree TBR (collapse / re-resolution).
#
# Discriminator: enumerate the FULL unrooted-TBR neighbourhood of the TS optimum
# with TreeSearch's SEPARATE, UNOPTIMISED R/Rcpp enumerator (TBR(tree, -1) from
# rearrange.cpp -- a different code path from ts_tbr.cpp).  If it finds an
# improving neighbour, the kernel missed a canonical TBR move => cause (i).
source("dev/benchmarks/tbr_shared_start_lib.R")

d <- prepareDataset("Zanol2014")
norm    <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))
asPhylo <- function(edge) structure(list(edge = edge, Nnode = d$nTip - 1L,
                          tip.label = names(d$phy)), class = "phylo")

RootInvariantTbr <- function(startTree, seed, rerootTips = names(d$phy)) {
  cur <- TsTbr(d, startTree, seed = seed, acceptEqual = FALSE)
  best <- cur$tree; bestLen <- cur$row$final_len
  repeat {
    improved <- FALSE
    for (tp in rerootTips) {
      rr <- norm(ape::root(best, outgroup = tp, resolve.root = TRUE))
      r  <- TsTbr(d, rr, seed = seed, acceptEqual = FALSE)
      if (r$row$final_len < bestLen) { bestLen <- r$row$final_len; best <- r$tree; improved <- TRUE }
    }
    if (!improved) break
  }
  list(len = bestLen, tree = norm(best))
}

# Reconstruct (and cache) the random-seed1 reroot-invariant optimum (~1284).
cacheFile <- "dev/benchmarks/tbr_results/ts_reroot_invariant_opt.tre"
if (file.exists(cacheFile)) {
  tsOpt <- norm(ape::read.tree(cacheFile)); tsLen <- TreeLength(tsOpt, d$phy)
  cat(sprintf("Loaded cached TS optimum: len=%.0f\n", tsLen))
} else {
  st <- norm({ set.seed(1001); RandomTree(d$phy, root = TRUE) })
  ri <- RootInvariantTbr(st, seed = 1)
  tsOpt <- ri$tree; tsLen <- ri$len
  ape::write.tree(tsOpt, cacheFile)
  cat(sprintf("Computed TS reroot-invariant optimum: len=%.0f (cached)\n", tsLen))
}

# Probe a tree's full unrooted-TBR neighbourhood with the R/Rcpp enumerator
# (TBRMoves -> all_tbr in rearrange.cpp -- a different code path from
# ts_tbr.cpp).  Scores every neighbour; reports the best + how many improve,
# and saves the best improving neighbour for move characterisation.
probeNeighbourhood <- function(tree, label, baseLen, saveBest = NULL) {
  cat(sprintf("\n== %s (len=%.0f) ==\n", label, baseLen))
  nb <- TBRMoves(norm(tree))
  n <- length(nb)
  ls <- vapply(nb, TreeLength, double(1), d$phy)
  best <- min(ls); bestIdx <- which.min(ls); nBetter <- sum(ls < baseLen - 0.5)
  cat(sprintf("  enumerated+scored %d TBR neighbours : best = %.0f  (%d improving)\n",
              n, best, nBetter))
  if (best < baseLen - 0.5) {
    cat(sprintf("  => KERNEL INCOMPLETE: canonical TBR improves it by %.0f (%.0f -> %.0f).\n",
                baseLen - best, baseLen, best))
    if (!is.null(saveBest)) {
      ape::write.tree(norm(nb[[bestIdx]]), saveBest)
      cat("     best improving neighbour saved ->", saveBest, "\n")
    }
  } else {
    cat("  => no single canonical-TBR move improves it (true TBR optimum by the enumerator).\n")
  }
  invisible(list(best = best, nBetter = nBetter, n = n))
}

# (1) The TS optimum the kernel converged on (all-tips). Does canonical TBR improve it?
probeNeighbourhood(tsOpt, "TS reroot-invariant optimum", tsLen,
                   saveBest = "dev/benchmarks/tbr_results/ts_opt_best_neighbour.tre")

# (2) Sanity: TNT's own optimum -- the enumerator should agree it's TBR-optimal
#     (TS kernel already holds it). Recompute a TNT optimum to probe.
tntRow <- TntTbr(d, norm({ set.seed(1001); RandomTree(d$phy, root = TRUE) }),
                 seed = 1, mulpars = FALSE, hold = 1)
tntOpt <- attr(tntRow, "tree")
if (!is.null(tntOpt))
  probeNeighbourhood(tntOpt, "TNT own optimum", tntRow$final_len)
