# Broad byte-identity gate for the incremental exact_verify rescore (TS_NA_INCR,
# commit 50046bfa + c78cc04d): merge-prep before default-on. For one bundled NA
# dataset (cell index), run a direct ts_tbr_search climb from several random
# starts in BOTH EW and IW regimes, with TS_NA_INCR off vs on, and assert the
# final score is byte-identical. Reports per-dataset pass count (or DIFF lines).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
                                              winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
nms <- names(inapplicable.phyData)
idx <- as.integer(if (length(commandArgs(trailingOnly = TRUE)))
                    commandArgs(trailingOnly = TRUE)[1]
                  else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
if (idx < 0L || idx >= length(nms)) stop("cell index out of range")
nm <- nms[idx + 1L]
phy <- inapplicable.phyData[[nm]]
ct <- attr(phy, "contrast"); lv <- attr(phy, "levels")
tip <- matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE)
wt <- attr(phy, "weight"); mins <- as.integer(MinimumLength(phy, compress = TRUE))
nm_lab <- names(phy)
starts <- c(1, 42, 100)
ntot <- 0L; nok <- 0L
for (conc in c(-1, 10)) {
  for (st in starts) {
    set.seed(st); redge <- RandomTree(nm_lab, root = TRUE)$edge
    run <- function(incr) {
      if (incr) Sys.setenv(TS_NA_INCR = "1") else Sys.unsetenv("TS_NA_INCR")
      set.seed(st)
      r <- TreeSearch:::ts_tbr_search(redge, ct, tip, wt, lv, maxHits = 1L,
             min_steps = mins, concavity = conc)
      Sys.unsetenv("TS_NA_INCR"); r$score
    }
    off <- run(FALSE); on <- run(TRUE)
    ntot <- ntot + 1L
    if (isTRUE(all.equal(off, on, tolerance = 0))) nok <- nok + 1L
    else cat(sprintf("DIFF %s conc=%s st=%d off=%.6f on=%.6f\n", nm, conc, st, off, on))
  }
}
cat(sprintf("IDENT %-16s %d/%d byte-identical (%d taxa)\n", nm, nok, ntot, length(phy)))
