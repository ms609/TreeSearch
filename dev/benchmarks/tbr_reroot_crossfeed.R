# tbr_reroot_crossfeed.R -- THE gating experiment (advisor, 2026-06-18).
#
# tbr_crossfeed.R only ever fed TNT the ROOTED 1302 optimum, which just
# re-proves root-dependence.  The decisive, never-run test is the
# REROOT-INVARIANT cross-feed:
#
#   (1) Take the all-tips reroot-invariant TS optimum (~1284 -- which SHOULD
#       be a complete unrooted-TBR local optimum, since every edge's large
#       side holds some tip so all-tip rooting overcomes the L812 skip).
#       Feed it into TNT bbreak.
#         TNT IMPROVES it  => residual is NEIGHBOURHOOD: TNT reaches moves our
#                             "complete" unrooted TBR doesn't (emulation not
#                             actually complete, or bbreak does more than
#                             single-tree TBR).  Fix THAT before building.
#         TNT HOLDS it     => residual is BASIN/PATH: 1284 and 1264 are both
#                             valid unrooted optima; TNT just navigates better.
#                             The lever is clip-ordering / restarts, NOT the
#                             move set.  Build still banks the root-dependence
#                             half, but won't reach TNT.
#   (2) Reciprocal: feed TNT's own ~1264 optimum into RootInvariantTbr.
#       Holds => confirms 1264 is a unrooted-TBR optimum too (basin/path).
source("dev/benchmarks/tbr_shared_start_lib.R")

d <- prepareDataset("Zanol2014")
norm    <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))
asPhylo <- function(edge) structure(list(edge = edge, Nnode = d$nTip - 1L,
                          tip.label = names(d$phy)), class = "phylo")

# All-tips reroot-invariant TBR to convergence.  Returns BOTH length and tree.
# (tbr_reroot_recovery.R's version returned only the length.)
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
  list(len = bestLen, tree = norm(best))
}

starts <- list(
  wagner = { set.seed(2001)
             w <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data,
                                                     d$weight, d$levels)
             norm(asPhylo(w$edge)) },
  random = norm({ set.seed(1001); RandomTree(d$phy, root = TRUE) }))

cat("=== Reroot-invariant cross-feed (Zanol2014, TNT target ~1261) ===\n\n")
for (sn in names(starts)) {
  st <- starts[[sn]]
  for (s in 1:2) {
    cat(sprintf("[%s seed=%d]  start_len=%.0f\n", sn, s, TreeLength(st, d$phy)))

    # (1) TS reroot-invariant optimum, then feed into TNT.
    ri <- RootInvariantTbr(st, seed = s, acceptEqual = FALSE)
    tnt_nomp <- TntTbr(d, ri$tree, seed = s, mulpars = FALSE, hold = 1)$final_len
    tnt_mp   <- TntTbr(d, ri$tree, seed = s, mulpars = TRUE,  hold = 1000)$final_len
    cat(sprintf("  TS reroot-invariant opt = %.0f\n", ri$len))
    cat(sprintf("    -> TNT bbreak nomulpars  : %.0f -> %.0f   (%s)\n",
                ri$len, tnt_nomp,
                if (tnt_nomp < ri$len - 0.5) "NEIGHBOURHOOD: TNT improves it"
                else "holds"))
    cat(sprintf("    -> TNT bbreak mulpars1000: %.0f -> %.0f\n", ri$len, tnt_mp))

    # (2) TNT's own optimum from the same start, fed back into reroot-invariant TS.
    tntRow <- TntTbr(d, st, seed = s, mulpars = FALSE, hold = 1)
    tntOpt <- attr(tntRow, "tree")
    if (!is.null(tntOpt)) {
      ri2 <- RootInvariantTbr(norm(tntOpt), seed = s, acceptEqual = FALSE)
      cat(sprintf("  TNT own opt = %.0f  -> reroot-invariant TS : %.0f -> %.0f   (%s)\n",
                  tntRow$final_len, tntRow$final_len, ri2$len,
                  if (ri2$len < tntRow$final_len - 0.5) "TS improves TNT"
                  else "TS HOLDS TNT (1264-class is a unrooted-TBR optimum)"))
    }
    cat("\n")
  }
}
cat("Reading: TNT improving the reroot-invariant opt => neighbourhood gap (fix first).\n",
    "TNT holding it while its own search still reaches lower => basin/path (clip order).\n",
    sep = "")
