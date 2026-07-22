# tbr_crossfeed.R -- the discriminating 2x2: feed each engine's local optimum
# into the other.  Distinguishes "TS reaches a worse basin" (escape/path
# problem) from "TS terminates before a true TBR local optimum"
# (neighbourhood incompleteness / premature stop).
source("dev/benchmarks/tbr_shared_start_lib.R")
d <- prepareDataset("Zanol2014")

set.seed(11)
wag <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
wagTree <- Preorder(RenumberTips(structure(list(edge=wag$edge, Nnode=d$nTip-1L,
                  tip.label=names(d$phy)), class="phylo"), names(d$phy)))
cat("Shared start = 1478\n\n")

# Each engine to its own local optimum from 1478.
tsRun  <- TsTbr(d, wagTree, seed = 2, acceptEqual = FALSE)
tsOpt  <- tsRun$tree
cat("TS strict descent  -> ", tsRun$row$final_len, " (TS local optimum)\n")

tntRow <- TntTbr(d, wagTree, seed = 2, mulpars = FALSE, hold = 1, randclip = TRUE)
tntOpt <- attr(tntRow, "tree")
cat("TNT bbreak nomulpars-> ", tntRow$final_len, " (TNT local optimum)\n\n")

# (1) TS local optimum -> TNT bbreak.  Does TNT improve a TS-converged tree?
cat("=== (1) Feed TS local optimum (", tsRun$row$final_len, ") into TNT bbreak ===\n", sep="")
for (rc in c(FALSE, TRUE)) {
  r <- TntTbr(d, tsOpt, seed = 2, mulpars = FALSE, hold = 1, randclip = rc)
  cat(sprintf("   TNT nomulpars randclip=%-5s : %s -> %s\n", rc, r$start_len, r$final_len))
}
rb <- TntTbr(d, tsOpt, seed = 2, mulpars = TRUE, hold = 1000, randclip = TRUE)
cat(sprintf("   TNT mulpars hold1000        : %s -> %s\n\n", rb$start_len, rb$final_len))

# (2) TNT local optimum -> TS strict bbreak.  Holds (escape) or wanders (bug)?
cat("=== (2) Feed TNT local optimum (", tntRow$final_len, ") into TS strict TBR ===\n", sep="")
r2 <- TsTbr(d, tntOpt, seed = 2, acceptEqual = FALSE)
cat(sprintf("   TS strict descent : %s -> %s  (converged=%s)\n",
            r2$row$start_len, r2$row$final_len, r2$converged))
r2e <- TsTbr(d, tntOpt, seed = 2, acceptEqual = TRUE, maxHits = 50L)
cat(sprintf("   TS plateau (mh50) : %s -> %s\n", r2e$row$start_len, r2e$row$final_len))
