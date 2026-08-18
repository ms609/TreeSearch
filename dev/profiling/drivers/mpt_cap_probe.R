# Decisive probe: can MaximizeParsimony return >256 distinct MPTs?
# Construct a soft polytomy: 4-taxon resolved backbone + an 8-taxon fan that
# is monophyletic but internally unresolved.  True resolved-MPT count =
# (2*8-3)!! = 135135  >>  256.  poolMaxSize = 2000.
#
# NOTE: with the 8-taxon fan this run is poolMaxSize-bound (returns ~2000 either
# way), so it does NOT, on its own, show the enum dedup-consistency fix or the
# collapse= option.  To see those, A/B the kill-switch on a SMALLER fan where the
# resolved count < poolMaxSize: TS_ENUM_RESOLVED=1 (old resolved dedup) vs unset
# (new collapsed dedup), e.g. fan=6 returns 947 vs 202; and compare
# collapse = FALSE vs collapse = TRUE (fan -> 1).  See memory
# mpt-enumeration-dedup-asymmetry.
suppressMessages(devtools::load_all(".", quiet = TRUE))

# ---- build matrix directly as a character matrix, then phyDat ----
backbone <- c("b1", "b2", "b3", "b4")
fan <- paste0("f", 1:8)
taxa <- c(backbone, fan)

mat <- matrix("0", nrow = length(taxa), ncol = 0,
              dimnames = list(taxa, NULL))

addchar <- function(m, ones) {
  col <- ifelse(rownames(m) %in% ones, "1", "0")
  cbind(m, col)
}
# Fan synapomorphy (x3 -> strongly supported fan clade)
for (i in 1:3) mat <- addchar(mat, fan)
# Backbone: (b1,b2) cherry x2, (b3,b4) cherry x2
for (i in 1:2) mat <- addchar(mat, c("b1", "b2"))
for (i in 1:2) mat <- addchar(mat, c("b3", "b4"))
# No character distinguishes f1..f8 -> internal fan is a soft polytomy.

phy <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))

set.seed(1L)
trees <- MaximizeParsimony(
  phy, concavity = Inf,
  maxReplicates = 20L, maxSeconds = 60,
  strategy = "intensive",
  control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 200L),
  verbosity = 2L)

cat("\n==== RESULT ====\n")
cat("returned trees:", length(trees), "\n")
cat("score attr:", attr(trees, "score"), "\n")
cat("n_topologies attr:", attr(trees, "n_topologies"), "\n")

# robust dedup on resolved unrooted topologies via Newick of sorted trees
canon <- vapply(trees, function(t) {
  ape::write.tree(TreeTools::SortTree(ape::unroot(t)))
}, character(1))
cat("distinct (sorted Newick):", length(unique(canon)), "\n")
