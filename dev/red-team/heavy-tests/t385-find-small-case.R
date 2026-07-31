# Find a SMALL xform dataset + tree whose Sankoff term genuinely varies with the
# rooting, for use as a self-guarding regression test.  Measured with
# ts_sankoff_test() directly -- the raw kernel, which the T-385 fix does not
# touch -- so the sensitivity is a property of the data, not of TreeLength().

lib <- if (dir.exists(".agent-t385")) ".agent-t385" else NULL
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = lib)
  library("TreeTools")
})

KernelSpread <- function(ds, h, tree) {
  recoded <- RecodeHierarchy(ds, h)
  xf <- TreeSearch:::.PrepareXformArgs(recoded, length(ds))
  vals <- vapply(names(ds), function(taxon) {
    tr <- RenumberTips(Renumber(RootTree(tree, taxon)), names(ds))
    TreeSearch:::ts_sankoff_test(tr[["edge"]], xf$n_states, xf$cost_matrices,
                                 xf$tip_states, xf$forced_root, xf$combo_grids,
                                 xf$tip_sec_known)$score
  }, numeric(1))
  vals
}

best <- NULL
for (seed in 1:400) {
  set.seed(seed)
  nTip <- 8L
  tipNames <- paste0("t", seq_len(nTip))
  primary <- sample(c("0", "1"), nTip, replace = TRUE)
  if (length(unique(primary)) < 2L) next
  sec1 <- ifelse(primary == "0", "-", sample(c("0", "1"), nTip, replace = TRUE))
  sec2 <- ifelse(primary == "0", "-", sample(c("0", "1"), nTip, replace = TRUE))
  mat <- cbind(primary, sec1, sec2)
  dimnames(mat) <- list(tipNames, NULL)
  ds <- try(phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                             ambiguity = "?"), silent = TRUE)
  if (inherits(ds, "try-error")) next
  h <- CharacterHierarchy("1" = 2:3)
  tree <- RandomTree(tipNames, root = TRUE)
  vals <- try(KernelSpread(ds, h, tree), silent = TRUE)
  if (inherits(vals, "try-error")) next
  sp <- diff(range(vals))
  if (sp > 0) {
    cat(sprintf("seed %3d: spread %g  values %s\n", seed, sp,
                paste(vals, collapse = " ")))
    if (is.null(best) || sp > best$spread) {
      best <- list(seed = seed, spread = sp, mat = mat, vals = vals,
                   newick = ape::write.tree(tree))
    }
    if (sp >= 2) break
  }
}

if (is.null(best)) {
  cat("no rooting-sensitive 8-tip case found\n")
} else {
  cat("\n=== BEST ===\nseed:", best$seed, " spread:", best$spread, "\n")
  cat("newick:", best$newick, "\n")
  cat("matrix:\n")
  print(best$mat)
  cat("\nvalues by root taxon:\n")
  print(best$vals)
}
