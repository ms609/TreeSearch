# Make-or-break: is the +356 TS Wagner deficit real, or a reconstruction
# artifact?  Compare the C++ kernel's OWN score to TreeLength of the
# reconstructed tree, and triangulate with the tested AdditionTree() path.
suppressMessages({
  library(TreeSearch, lib.loc = "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem")
  library(TreeTools)
})
nm   <- Sys.getenv("DS", "Zanol2014")
phy  <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
taxa <- names(phy); n <- length(taxa)
at <- attributes(phy)
contrast <- at$contrast
tipData  <- matrix(unlist(phy, use.names = FALSE), nrow = n, byrow = TRUE)
weight   <- TreeSearch:::.ScaleWeight(at$weight)
levels   <- at$levels

.EdgeToPhylo <- function(edge) {
  tr <- structure(list(edge = edge, tip.label = taxa, Nnode = n - 1L), class = "phylo")
  Renumber(tr)
}

cat("== kernel score vs TreeLength(reconstruction) ==\n")
for (i in 1:5) {
  set.seed(i)
  res <- TreeSearch:::ts_random_wagner_tree(contrast, tipData, weight, levels)
  tr  <- .EdgeToPhylo(res$edge)
  tl  <- TreeLength(tr, phy)
  cat(sprintf("  seed %d: kernel=%.0f  TreeLength=%.0f  %s\n",
              i, res$score, tl,
              if (abs(res$score - tl) < 0.5) "MATCH" else "*** MISMATCH ***"))
}

cat("\n== AdditionTree() (tested path) with random sequences ==\n")
for (i in 1:5) {
  set.seed(100 + i)
  seq_i <- sample(taxa)
  tr <- AdditionTree(phy, sequence = seq_i)
  cat(sprintf("  seed %d: AdditionTree TreeLength=%.0f\n", i, TreeLength(tr, phy)))
}

cat("\n== sanity: a purely random topology score (upper reference) ==\n")
set.seed(7)
rt <- RandomTree(phy, root = taxa[1])
cat(sprintf("  RandomTree TreeLength=%.0f\n", TreeLength(rt, phy)))
