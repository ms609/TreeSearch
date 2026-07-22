# tbr_unrooted_validate.R -- quality + PERF of the in-kernel unrooted TBR
# (TBRParams::unrooted, reroot-at-convergence) on real data (Zanol2014).
#
# Correctness (0 canonical-improving) is proven separately by the small-tree
# oracle (tbr_oracle.R, unrooted=1 emul=0).  This script measures, per start:
#   - final length   rooted (default) vs unrooted (in-kernel reroot)
#   - wall-clock time rooted vs unrooted  => the per-tbr_search perf cost
# from both POOR (random) and GOOD (RAS-Wagner) starts.  Context: TNT reaches
# ~1262-1264; closing to TRUE unrooted-TBR optima (this fix) is expected to land
# ~1265-1272 -- the residual to TNT is basin/escape, a SEPARATE mechanism.
source("dev/benchmarks/tbr_shared_start_lib.R")
d <- prepareDataset("Zanol2014")
norm <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))
asPhylo <- function(edge) structure(list(edge = edge, Nnode = d$nTip - 1L,
                          tip.label = names(d$phy)), class = "phylo")

runKernel <- function(tree, seed, unrooted) {
  edge <- PhyloToKernelEdge(tree, d)
  set.seed(seed)
  t <- system.time(
    res <- TreeSearch:::ts_tbr_diagnostics(
      edge, d$contrast, d$tip_data, d$weight, d$levels,
      maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, unrooted = unrooted))
  tr <- norm(asPhylo(res$edge))
  list(tree = tr, len = TreeLength(tr, d$phy), sec = as.double(t["elapsed"]))
}

# Is `tree` canonical-unrooted-TBR clean? all_tbr at two rootings (covers all
# break edges).  Expensive (~2x100k neighbours); call sparingly.
isClean <- function(tree) {
  base <- TreeLength(tree, d$phy)
  best <- base
  for (rt in names(d$phy)[1:2]) {
    nb <- TBRMoves(norm(RootTree(tree, rt)))
    best <- min(best, min(vapply(nb, TreeLength, double(1), d$phy)))
  }
  best >= base - 0.5
}

mkStart <- function(kind, seed) {
  if (kind == "random") norm({ set.seed(1000 + seed); RandomTree(d$phy, root = TRUE) })
  else { set.seed(2000 + seed)
         w <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
         norm(asPhylo(w$edge)) }
}

cat("=== In-kernel unrooted TBR: quality + perf (Zanol2014, n=74; TNT ~1262-1264) ===\n")
cat(sprintf("%-7s %-4s %-7s | %-8s %-7s | %-8s %-7s %-6s | %-5s\n",
            "start","seed","startL","rootedL","sec","unrootL","sec","clean","x"))
rows <- list()
for (kind in c("random","wagner")) for (s in 1:3) {
  st <- mkStart(kind, s); sl <- TreeLength(st, d$phy)
  r <- runKernel(st, s, FALSE)
  u <- runKernel(st, s, TRUE)
  # Cleanliness on full 74-tip is ~330s/call; correctness is proven broadly by
  # the small-tree oracle, so confirm on real 74-tip data for ONE start only.
  clean <- if (kind == "random" && s == 1) isClean(u$tree) else NA
  cat(sprintf("%-7s %-4d %-7.0f | %-8.0f %-7.2f | %-8.0f %-7.2f %-6s | %-5.1f\n",
              kind, s, sl, r$len, r$sec, u$len, u$sec, clean, u$sec / r$sec))
  rows[[length(rows)+1]] <- data.frame(kind, seed=s, startL=sl,
      rootedL=r$len, rootedSec=r$sec, unrootL=u$len, unrootSec=u$sec,
      clean=clean, ratio=u$sec/r$sec)
}
res <- do.call(rbind, rows)
cat(sprintf("\nMEDIAN: rooted=%.0f  unrooted=%.0f  (gain %.0f)   median time x%.1f\n",
            median(res$rootedL), median(res$unrootL),
            median(res$rootedL - res$unrootL), median(res$ratio)))
cat(sprintf("unrooted results canonical-TBR-clean: %d/%d\n", sum(res$clean), nrow(res)))
write.csv(res, "dev/benchmarks/tbr_results/tbr_unrooted_validate.csv", row.names = FALSE)
