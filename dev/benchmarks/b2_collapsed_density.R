# B2 collapsed-density probe (2026-06-22, architecture-audit B2 follow-up).
#
# CONTEXT: the full-polytomy plan (2026-03-22) proposed building collapsed-region
# regraft merging as "the main win" (Phase 3, 3-5 days).  CODE READ shows it is
# ALREADY IMPLEMENTED AND ACTIVE in production hill-climbing / sectorial / ratchet
# TBR (ts_tbr.cpp:1581 `if (use_collapsed && collapsed[below]) continue;` +
# :1692 reroot kept_ei filter; use_collapsed is ON whenever collect_pool==nullptr,
# which is every hill-climbing pass -- only the final MPT-enumeration pass disables
# it).  So Phase 3 is not a build opportunity; it is running.
#
# The only open question for the WALL-CLOCK mission ([[mission-wallclock-to-optimum]]):
# does that active merging actually skip a meaningful fraction of work, or is the
# morphological data class too barren of zero-length edges (the plan's own "0%
# skip rate" note) for it to matter?  Both clip-skip (Phase 2, counted by
# n_zero_skipped) and regraft-merge (Phase 3, uncounted) key on the SAME collapsed[]
# array, so n_zero_skipped/n_evaluated over a real descent is a faithful proxy for
# the collapsed-edge density where wall-clock is spent.
#
# This is a NO-BUILD measurement: ts_tbr_diagnostics already returns n_zero_skipped
# and n_evaluated, accumulated across all passes of the descent (work-weighted).
# We run a full TBR descent from a POOR start (random tree = upper bound on
# collapsed density -- a random start has more unsupported edges than a Wagner
# start, so this is the case most favourable to B2) on Zanol2014.
#
# Read: if the skip rate is ~0% even from a random start, the data class is barren
# of zero-length edges => the active merging captures ~nothing and there is no
# wall-clock win in B2's main mechanism for this class.  If it is substantial,
# the win is already being CAPTURED (merging is on) => still nothing new to build,
# but worth confirming the trajectory benefit is realised.
#
# Usage: TS_LIB=.agent-fuse Rscript dev/benchmarks/b2_collapsed_density.R

suppressMessages({
  library(TreeSearch,
          lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-fuse"),
                                  winslash = "/", mustWork = TRUE))
  library(TreeTools)
})

.FitchPhy <- function(p) {
  m <- PhyDatToMatrix(p, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}
e <- new.env()
utils::data("inapplicable.phyData", package = "TreeSearch", envir = e)
phy <- .FitchPhy(e[["inapplicable.phyData"]][["Zanol2014"]])
at <- attributes(phy)
labels <- names(phy)
nTip <- length(phy)
tipData <- matrix(unlist(phy, use.names = FALSE), nrow = nTip, byrow = TRUE)

cat(sprintf("Zanol2014: %d tips, %d patterns\n", nTip, length(at$weight)))

# Run a full in-kernel TBR descent from several poor (random) starts and read the
# accumulated diagnostics.  Random start = high collapsed density (upper bound).
runOne <- function(seed) {
  set.seed(seed)
  start <- RandomTree(phy, root = TRUE)
  edge <- Preorder(RenumberTips(start, labels))[["edge"]]
  set.seed(seed)
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, at$contrast, tipData, at$weight, at$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L)
  res
}

r1 <- runOne(1L)
cat("\nReturned diagnostic fields:\n")
str(r1[setdiff(names(r1), "edge")])

# Aggregate over a few starts.
seeds <- 1:5
rs <- lapply(seeds, runOne)
getf <- function(r, f) if (!is.null(r[[f]])) as.double(r[[f]]) else NA_real_
zskip <- vapply(rs, getf, double(1), "n_zero_skipped")
neval <- vapply(rs, getf, double(1), "n_evaluated")

cat("\n=== Collapsed clip-skip over full TBR descents (random starts) ===\n")
for (i in seq_along(seeds)) {
  cat(sprintf("  seed %d: n_zero_skipped=%g  n_evaluated=%g  skipRate=%.3f%%\n",
              seeds[i], zskip[i], neval[i],
              100 * zskip[i] / max(neval[i] + zskip[i], 1)))
}
totalSkip <- sum(zskip, na.rm = TRUE)
totalEval <- sum(neval, na.rm = TRUE)
cat(sprintf("\nPooled clip-skip rate: %g / %g = %.3f%% of clip-candidate work\n",
            totalSkip, totalEval + totalSkip,
            100 * totalSkip / max(totalEval + totalSkip, 1)))
cat("(Phase-2 clip-skip; Phase-3 regraft-merge keys on the same collapsed[] array,\n",
    " so this proxies the collapsed-edge density where wall-clock is spent.)\n", sep = "")
