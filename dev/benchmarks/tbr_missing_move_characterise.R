# tbr_missing_move_characterise.R -- runs only if the neighbourhood probe shows
# the kernel is incomplete (canonical TBR improves the TS reroot-invariant
# optimum).  Pins WHICH kernel pruning drops the improving move, no rebuild.
#
# Suspects in ts_tbr.cpp:
#   * L812 smaller-subtree skip  -> improving move clips the LARGER side only.
#   * L817/L919 collapsed pruning -> improving move touches a ZERO-LENGTH edge.
#   * indirect-scoring cutoff / vp-dedup -> neither of the above.
source("dev/benchmarks/tbr_shared_start_lib.R")
d <- prepareDataset("Zanol2014")
norm <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))

optFile <- "dev/benchmarks/tbr_results/ts_reroot_invariant_opt.tre"
nbFile  <- "dev/benchmarks/tbr_results/ts_opt_best_neighbour.tre"
stopifnot(file.exists(optFile), file.exists(nbFile))
opt <- norm(ape::read.tree(optFile))
nb  <- norm(ape::read.tree(nbFile))
optLen <- TreeLength(opt, d$phy); nbLen <- TreeLength(nb, d$phy)
cat(sprintf("TS optimum = %.0f ; best canonical-TBR neighbour = %.0f (improve %.0f)\n\n",
            optLen, nbLen, optLen - nbLen))

# --- Zero-length edges of the optimum (the collapsed-pruning suspects) ---
# An internal edge is "zero length" if collapsing it (merging child into parent)
# leaves the parsimony length unchanged.
ed <- opt[["edge"]]; nTip <- d$nTip
internalChildEdges <- which(ed[, 2] > nTip)   # edges whose child is internal
zeroLen <- 0L
for (e in internalChildEdges) {
  collapsed <- opt
  collapsed$edge.length <- NULL
  collapsed <- ape::di2multi(  # collapse just this edge via a tiny length vector
    { t2 <- opt; t2$edge.length <- rep(1, nrow(ed)); t2$edge.length[e] <- 0; t2 },
    tol = 0.5)
  if (abs(TreeLength(collapsed, d$phy) - optLen) < 0.5) zeroLen <- zeroLen + 1L
}
cat(sprintf("zero-length internal edges in the optimum: %d / %d\n",
            zeroLen, length(internalChildEdges)))
cat(if (zeroLen == 0)
      "  => collapsed-edge pruning is INACTIVE here; cause is cutoff/dedup, not collapsed.\n"
    else
      "  => collapsed-edge pruning is a LIVE suspect (zero-length edges present).\n")

# --- Move magnitude: splits that differ between optimum and best neighbour ---
spOpt <- TreeTools::as.Splits(opt)
spNb  <- TreeTools::as.Splits(nb,  tipLabels = TipLabels(opt))
# Count splits in one but not the other (RF-style raw difference).
inOpt <- apply(as.logical(spOpt), 1, function(r) paste(as.integer(r), collapse = ""))
inNb  <- apply(as.logical(spNb),  1, function(r) paste(as.integer(r), collapse = ""))
# Normalise complement (a split and its complement are the same bipartition).
canon <- function(s) { v <- as.integer(strsplit(s, "")[[1]])
  if (v[1] == 1) paste(1L - v, collapse = "") else s }
inOptC <- vapply(inOpt, canon, ""); inNbC <- vapply(inNb, canon, "")
nDiff <- length(setdiff(inOptC, inNbC))
cat(sprintf("\nsplits in optimum absent from best neighbour: %d (TBR move magnitude)\n", nDiff))

cat("\nReading:\n",
    " - zero-length edges present + small splits-diff => collapsed pruning the likely culprit.\n",
    " - no zero-length edges => indirect-cutoff / vp-dedup drops the move.\n", sep = "")
