# Is our default TBR's "convergence" genuine, on the SHIPPING cpp-search build
# (post directional-vroot fix, commit 2b299e4b)?  We run TBR to convergence via
# ts_tbr_diagnostics, then enumerate the FULL unrooted canonical-TBR neighbourhood
# of the result with the SEPARATE, unoptimised enumerator TBRMoves (-> all_tbr in
# rearrange.cpp, a different code path).  >0 improving neighbour => the kernel
# falsely declared convergence (the competent-chaum move-incompleteness finding).
#
# Result (Zanol2014, 2026-06-18): GOOD Wagner starts -> genuine optima (0
# improving); POOR random starts -> strand at 1272 with only 1-9 improving (vs the
# chip's PRE-fix 40+), i.e. the vroot scoring fix recovered most of the gap and the
# residual move-skip bug is small.  See dev/plans/2026-06-18-wagner-insertion-cost-bug.md.
#
# Env: TS_LIB (default .agent-wagsect).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-wagsect"),
            winslash = "/"))
  library(TreeTools)
})

data("inapplicable.phyData", package = "TreeSearch")
fitchify <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy <- fitchify(inapplicable.phyData[["Zanol2014"]])
at  <- attributes(phy)
d <- list(phy = phy, contrast = at$contrast,
          tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
          weight = at$weight, levels = at$levels, nTip = length(phy))
norm <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))

# Run our default (rooted, optimised) TBR to convergence from a warm start.
tsTbr <- function(start, seed) {
  set.seed(seed)
  res <- TreeSearch:::ts_tbr_diagnostics(norm(start)[["edge"]], d$contrast, d$tip_data,
           d$weight, d$levels, maxHits = 1L, acceptEqual = FALSE)
  structure(list(edge = res$edge, Nnode = d$nTip - 1L, tip.label = names(d$phy)),
            class = "phylo")
}

probe <- function(label, start, seed) {
  tr <- tsTbr(start, seed)
  baseLen <- TreeLength(tr, d$phy)
  ls <- vapply(TBRMoves(norm(tr)), TreeLength, double(1), d$phy)   # full unrooted-TBR neighbourhood
  cat(sprintf("%-12s start=%4.0f | TS-TBR converged=%.0f | enum %d nb, best=%.0f, %d IMPROVING\n",
              label, TreeLength(norm(start), d$phy), baseLen, length(ls), min(ls),
              sum(ls < baseLen - 0.5)))
}

cat("--- GOOD starts (RAS Wagner, post-fix) ---\n")
for (s in 1:3) { set.seed(s); probe(sprintf("wagner s%d", s),
  AdditionTree(d$phy, sequence = sample(seq_along(d$phy))), s) }

cat("--- POOR starts (random topology, like the chip's) ---\n")
for (s in 1:3) { set.seed(1000 + s); probe(sprintf("random s%d", s),
  RandomTree(d$phy, root = TRUE), s) }
