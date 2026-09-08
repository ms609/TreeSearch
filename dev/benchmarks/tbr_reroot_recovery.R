# tbr_reroot_recovery.R -- Step-2 confirmation + calibration.
#
# If TS's deficit is a root-dependent (rooted) TBR neighbourhood, then emulating
# root-invariance with the EXISTING rooted kernel (TBR -> reroot-sweep -> TBR,
# looped to convergence over ALL tips) should recover a large part of the gap to
# TNT.  It does -- ~half -- but STALLS +15-36 above TNT, and plateau-crossing
# does not close the residual.  So root-dependence is a proven MAJOR contributor,
# not (yet) shown to be the whole cause.  See dev/plans/2026-06-18-tbr-shared-start.md.
source("dev/benchmarks/tbr_shared_start_lib.R")

d <- prepareDataset("Zanol2014")
norm    <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))
asPhylo <- function(edge) structure(list(edge = edge, Nnode = d$nTip - 1L,
                          tip.label = names(d$phy)), class = "phylo")

# Root-invariant TBR emulation: run rooted TBR; then try re-rooting at each tip
# in `rerootTips` and re-running; adopt any rooting that improves; repeat until a
# full reroot sweep yields no improvement.  Uses only the shipping rooted kernel.
RootInvariantTbr <- function(startTree, seed, acceptEqual = FALSE, maxHits = 1L,
                             rerootTips = names(d$phy)) {
  cur <- TsTbr(d, startTree, seed = seed, acceptEqual = acceptEqual, maxHits = maxHits)
  best <- cur$tree; bestLen <- cur$row$final_len
  repeat {
    improved <- FALSE
    for (tp in rerootTips) {
      rr <- norm(ape::root(best, outgroup = tp, resolve.root = TRUE))
      r  <- TsTbr(d, rr, seed = seed, acceptEqual = acceptEqual, maxHits = maxHits)
      if (r$row$final_len < bestLen) {
        bestLen <- r$row$final_len; best <- r$tree; improved <- TRUE
      }
    }
    if (!improved) break
  }
  bestLen
}

set.seed(2001)
w <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
starts <- list(wagner = norm(asPhylo(w$edge)),
               random = norm({set.seed(1001); RandomTree(d$phy, root = TRUE)}))

cat("Full 74-tip reroot-invariant TBR to convergence (Zanol):\n")
cat(sprintf("%-7s %-5s  %-10s %-12s %-12s %-6s\n",
            "start", "seed", "TS_rooted", "reroot_strict", "reroot_plateau", "TNT"))
for (sn in names(starts)) {
  st <- starts[[sn]]
  for (s in 1:2) {
    rooted  <- TsTbr(d, st, seed = s, acceptEqual = FALSE)$row$final_len
    riStr   <- RootInvariantTbr(st, seed = s, acceptEqual = FALSE)
    riPlat  <- RootInvariantTbr(st, seed = s, acceptEqual = TRUE, maxHits = 50L)
    tnt     <- TntTbr(d, st, seed = s, mulpars = FALSE, hold = 1)$final_len
    cat(sprintf("%-7s %-5d  %-10.0f %-12.0f %-12.0f %-6.0f\n",
                sn, s, rooted, riStr, riPlat, tnt))
  }
}
cat("\nReading: reroot-invariance recovers ~half the gap but stalls +15-36 above\n",
    "TNT; plateau-crossing does not close the residual.\n", sep = "")
