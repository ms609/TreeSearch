# tbr_pilot.R -- pre-flight gates #1 (length identity) + #2 (seed sensitivity,
# first/best-improving characterization) on ONE Zanol start tree, both modes.
source("dev/benchmarks/tbr_shared_start_lib.R")

d  <- prepareDataset("Zanol2014")

# A deliberately POOR start so TBR has room to climb.
set.seed(11)
wag <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
wagTree <- structure(list(edge = wag$edge, Nnode = d$nTip - 1L,
                          tip.label = names(d$phy)), class = "phylo")
wagTree <- Preorder(RenumberTips(wagTree, names(d$phy)))
cat("Start TreeLength(Wagner seed11) =", TreeLength(wagTree, d$phy), "\n\n")

seeds <- c(1, 2, 3)

cat("=== MODE A: strict descent (TS acceptEqual=F ; TNT nomulpars hold 1) ===\n")
rows <- list()
for (s in seeds) {
  rows[[length(rows)+1]] <- TntTbr(d, wagTree, seed = s, mulpars = FALSE, hold = 1)
  rows[[length(rows)+1]] <- TsTbr(d, wagTree, seed = s, acceptEqual = FALSE)$row
}
modeA <- do.call(rbind, rows)
print(modeA, row.names = FALSE)

cat("\n=== MODE B: buffer/plateau (TS acceptEqual=T ; TNT mulpars hold 1000) ===\n")
rows <- list()
for (s in seeds) {
  rows[[length(rows)+1]] <- TntTbr(d, wagTree, seed = s, mulpars = TRUE, hold = 1000)
  rows[[length(rows)+1]] <- TsTbr(d, wagTree, seed = s, acceptEqual = TRUE,
                                  maxHits = 5L)$row
}
modeB <- do.call(rbind, rows)
print(modeB, row.names = FALSE)

cat("\n--- GATE CHECKS ---\n")
cat("TNT start_len (R) vs start_len (TNT stdout):\n")
print(unique(modeA[, c("start_len", "start_len_tnt")]))
cat("TS final_len (R) vs final_len (kernel res$score) should match exactly:\n")
tsRows <- rbind(modeA[modeA$engine=="TS",], modeB[modeB$engine=="TS",])
print(tsRows[, c("seed","final_len","final_len_tnt")])
