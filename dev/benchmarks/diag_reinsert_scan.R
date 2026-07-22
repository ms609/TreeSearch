# Per-edge exactness probe: clip one tip from a fixed tree, and for every
# reattachment edge compare the union-of-finals formula cost and the directional
# intersect-else-union cost against the TRUE full-rescore cost.
suppressMessages({
  library(TreeSearch, lib.loc = "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem")
  library(TreeTools)
})
nm   <- Sys.getenv("DS", "Zanol2014")
phy  <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
n    <- length(phy)

set.seed(11)
tr <- AdditionTree(phy)                       # any valid full tree
phyO <- phy[tr$tip.label]                      # data in the tree's tip order
at <- attributes(phyO)
contrast <- at$contrast
tipData  <- matrix(unlist(phyO, use.names = FALSE), nrow = n, byrow = TRUE)
weight   <- TreeSearch:::.ScaleWeight(at$weight); levels <- at$levels

# --- small partial tree (12 taxa): is dir exact when sets are ambiguous? ---
{
  small <- names(phy)[1:12]
  phyS <- phy[small]; nS <- length(small)
  trS <- AdditionTree(phyS)
  phySO <- phyS[trS$tip.label]; atS <- attributes(phySO)
  tdS <- matrix(unlist(phySO, use.names = FALSE), nrow = nS, byrow = TRUE)
  rS <- TreeSearch:::ts_reinsert_scan(trS$edge, atS$contrast, tdS,
                                      TreeSearch:::.ScaleWeight(atS$weight), atS$levels, 4L)
  cat(sprintf("-- SMALL 12-taxon tree, clip tip 4: union==actual %d/%d, dir==actual %d/%d --\n",
              sum(rS$union_extra == rS$actual_extra), length(rS$actual_extra),
              sum(rS$dir_extra == rS$actual_extra), length(rS$actual_extra)))
}

for (clip in c(3L, 12L, 40L)) {
  r <- TreeSearch:::ts_reinsert_scan(tr$edge, contrast, tipData, weight, levels, clip)
  un <- r$union_extra; di <- r$dir_extra; ac <- r$actual_extra
  cat(sprintf("\n-- clip tip %d (%s), main_score=%d, %d edges --\n",
              clip, tr$tip.label[clip], r$main_score, length(ac)))
  cat(sprintf("  union==actual: %d/%d   | dir==actual: %d/%d\n",
              sum(un == ac), length(ac), sum(di == ac), length(ac)))
  cat(sprintf("  min actual=%d | union picks edge w/ actual=%d | dir picks edge w/ actual=%d\n",
              min(ac), ac[which.min(un)], ac[which.min(di)]))
  # show a few rows where they disagree with truth
  bad <- which(un != ac | di != ac)
  if (length(bad)) {
    show <- head(bad, 6)
    for (i in show) cat(sprintf("    edge(%d,%d): union=%d dir=%d actual=%d\n",
                                r$above[i], r$below[i], un[i], di[i], ac[i]))
  }
}
