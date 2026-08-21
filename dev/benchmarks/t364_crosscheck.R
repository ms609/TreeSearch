#!/usr/bin/env Rscript
# Cross-check: does THIS harness reproduce the numbers already on the record?
#
# Two recorded figures, both on cases outside the MorphoBank corpus, so agreement
# ties the battery's machinery to the existing findings rather than asking the
# reader to trust a new pipeline:
#   * T-364: the 6-tip letters dataset with constraint {e,f} violates in
#     35/400 = 8.75% of seeded addition orders, pre-fix.
#   * T-384: on a 22-tip matrix with a 6-taxon constraint, 19/200 = 9.5% of
#     arm-2 starts display the split only through the complement.
#
# The seed loop mirrors the recorded protocol exactly (set.seed(i) then let
# `sequence` default to random), rather than supplying explicit orders, because
# the recorded rate is a property of that protocol.
#
# Env: TS_LIB, ARM, COMMON, N_SEED (default 400), OUT_DIR.

suppressMessages({
  ts_lib <- Sys.getenv("TS_LIB", "")
  if (nzchar(ts_lib)) {
    library(TreeSearch, lib.loc = normalizePath(ts_lib, winslash = "/", mustWork = TRUE))
  } else {
    library(TreeSearch)
  }
  library(TreeTools)
})
source(Sys.getenv("COMMON", "t364_common.R"))

arm <- Sys.getenv("ARM", "?")
n_seed <- as.integer(Sys.getenv("N_SEED", "400"))
out_dir <- Sys.getenv("OUT_DIR", "t364_xcheck")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sweep <- function(label, dataset, cons_pd, group, n) {
  tip0 <- names(dataset)[1]
  cc <- cx <- ok <- 0L
  for (i in seq_len(n)) {
    set.seed(i)
    tr <- suppressWarnings(AdditionTree(dataset, constraint = cons_pd))
    cl <- classify_tree(tr, group, tip0)
    cc <- cc + cl$clade_canonical
    cx <- cx + cl$clade_complement
    ok <- ok + cl$compliant
  }
  viol <- n - ok
  comp_only <- ok - cc
  cat(sprintf("  %-22s arm %s: violating=%d/%d (%.2f%%)  complementONLY=%d/%d (%.2f%%)\n",
              label, arm, viol, n, 100 * viol / n, comp_only, n, 100 * comp_only / n))
  data.frame(arm = arm, case = label, n = n, n_violating = viol,
             pct_violating = 100 * viol / n, n_complement_only = comp_only,
             pct_complement_only = 100 * comp_only / n, stringsAsFactors = FALSE)
}

rows <- list()

# ---- Case 1: the T-364 test case, verbatim (test-AdditionTree.R:66-78).
ds6 <- TreeTools::MatrixToPhyDat(matrix(
  c(0, 1, 1, 1, 0, 1,
    0, 1, 1, 0, 0, 1), ncol = 2,
  dimnames = list(letters[1:6], NULL)))
cons6 <- c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1)
rows[[1]] <- sweep("T-364 6-tip {e,f}", ds6,
                   TreeTools::MatrixToPhyDat(cons6), c("e", "f"), n_seed)

# ---- Case 2: T-384's vehicle -- 22-tip congreveLamsdellMatrices[[1]], a 6-taxon
# ---- constraint.  The recorded row does not pin WHICH six taxa, so this is "a"
# ---- 6-taxon constraint (the first six tips), not a claim of verbatim identity.
data("congreveLamsdellMatrices", package = "TreeSearch")
ds22 <- congreveLamsdellMatrices[[1]]
tips22 <- names(ds22)
grp22 <- tips22[1:6]
cat(sprintf("  [case 2] 22-tip, group = %s (contains tip0 '%s': %s)\n",
            paste(grp22, collapse = ","), tips22[1], tips22[1] %in% grp22))
rows[[2]] <- sweep("T-384 22-tip first6", ds22,
                   constraint_phydat(grp22, tips22), grp22, min(n_seed, 200L))

# The first-six group CONTAINS tip 0, so canonicalisation makes the canonical mask
# the 16-tip side -- the `in0` geometry, where the complement rate is high.  Also
# run a group that EXCLUDES tip 0, which is the `no0` geometry the recorded 9.5%
# came from.
grp22b <- tips22[2:7]
cat(sprintf("  [case 2b] 22-tip, group = %s (contains tip0: %s)\n",
            paste(grp22b, collapse = ","), tips22[1] %in% grp22b))
rows[[3]] <- sweep("T-384 22-tip tips2-7", ds22,
                   constraint_phydat(grp22b, tips22), grp22b, min(n_seed, 200L))

D <- do.call(rbind, rows)
write.csv(D, file.path(out_dir, sprintf("xcheck_arm%s.csv", arm)), row.names = FALSE)
