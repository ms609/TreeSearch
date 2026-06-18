# tbr_verify.R -- sanity checks on the striking pilot result.
source("dev/benchmarks/tbr_shared_start_lib.R")
d <- prepareDataset("Zanol2014")
set.seed(11)
wag <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
wagTree <- Preorder(RenumberTips(structure(list(edge = wag$edge, Nnode = d$nTip-1L,
                  tip.label = names(d$phy)), class = "phylo"), names(d$phy)))
cat("start =", TreeLength(wagTree, d$phy), "\n\n")

# (1) Does TS TBR converge, and what does the per-pass trajectory look like?
for (ae in c(FALSE, TRUE)) {
  r <- TsTbr(d, wagTree, seed = 2, acceptEqual = ae, maxHits = if (ae) 5L else 1L)
  p <- r$passes
  cat(sprintf("TS acceptEqual=%-5s  final=%.0f  converged=%s  n_accepted=%d  n_passes=%d\n",
              ae, r$row$final_len, r$converged, r$n_accepted, nrow(p)))
  cat("  productive passes:", sum(p$productive), " null passes:", sum(!p$productive), "\n")
}

# (2) TNT determinism: norandclip same seed twice (should be identical);
#     randclip different seeds (should differ).
cat("\n--- TNT norandclip x2 (determinism) ---\n")
a1 <- TntTbr(d, wagTree, seed=1, mulpars=FALSE, hold=1, randclip=FALSE)
a2 <- TntTbr(d, wagTree, seed=1, mulpars=FALSE, hold=1, randclip=FALSE)
cat("norandclip seed1 run1:", a1$final_len, " run2:", a2$final_len, "\n")
b1 <- TntTbr(d, wagTree, seed=1, mulpars=FALSE, hold=1, randclip=TRUE)
b2 <- TntTbr(d, wagTree, seed=2, mulpars=FALSE, hold=1, randclip=TRUE)
cat("randclip seed1:", b1$final_len, " seed2:", b2$final_len, "\n")

# (3) Sanity: TNT bbreak from the OPTIMAL T0 (1271) must NOT do RAS (stay <=1271).
cat("\n--- TNT bbreak from T0=1271 (must not re-randomise) ---\n")
t0 <- ape::read.tree(file.path(T0_DIR, "Zanol2014.tre"))
c1 <- TntTbr(d, t0, seed=1, mulpars=FALSE, hold=1, randclip=TRUE)
cat("T0 start:", c1$start_len, " final:", c1$final_len, "\n")

# (4) TS from T0=1271 strict descent (should stay near 1271).
cat("\n--- TS bbreak from T0=1271 ---\n")
t0r <- TsTbr(d, t0, seed=1, acceptEqual=FALSE)
cat("T0 start: 1271  TS final:", t0r$row$final_len, "\n")
