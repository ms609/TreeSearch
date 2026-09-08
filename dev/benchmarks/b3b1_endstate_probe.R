# B3xB1 FALSIFICATION PROBE (2026-06-22, architecture-audit follow-up).
#
# Question: our highest-value mode, SECTORIAL search, runs TBR with the
# unrooted-completeness mechanism GATED OFF (ts_tbr.h:39-44 / ts_tbr.cpp:1384:
# do_reroot is false whenever sector_mask/cd/tabu/pool is active, because the
# outer reroot loop permutes node ids and would invalidate the node-id-keyed
# sector_mask).  So sectorial TBR can stop at a rooted-only local optimum that a
# root-edge / reroot move would escape.  Does that cost FINAL quality on the one
# open quality item -- reliably reaching 1261 on Zanol2014?
#
# Decisive number (advisor-framed, END-STATE not mid-sectorial): do the final
# converged Zanol replicates that MISS 1261 have an improving unrooted-TBR
# neighbour?
#   * 0-improving -> those misses are already TRUE unrooted-TBR optima (just the
#     wrong basin).  Completeness cannot rescue a locally-complete tree.  The
#     B3xB1 mechanism ("sectorial stops us ABOVE the unrooted optimum") is
#     refuted -> DEAD for the mission; pivot to B2 or close the frontier.
#   * >0-improving -> completeness is leaving quality on the table at end-state
#     -> earn the attribution work (which phase; trailing-TBR gate; sector remap).
#
# Why end-state, not mid-sectorial: sectorial output is mid-pipeline; ratchet/
# drift/fuse and trailing TBRs run after it and a plain-mode (gate-ON) TBR could
# recover whatever sectorial left.  The end-state count is the one that maps to
# the mission.  Prior is LOW: (a) main-path completeness "does NOT close the TNT
# level gap"; (b) good-Wagner-start production "already reaches 0-improving";
# (c) B1 attributed the reliable-1261 residual to budget, not a ceiling.  Run as
# falsification.
#
# Primary test = the kernel's OWN validated unrooted path
# (ts_tbr_diagnostics(unrooted=TRUE), the all-tips path that scored 0/60 on the
# completeness oracle) applied to each missed end-state tree.  TBRMoves() is
# infeasible at 74 tips (O(n^3) neighbourhood); SPRMoves() (O(n^2)) is a cheap
# cross-check that any kernel improvement is a genuine canonical move.
#
# Usage: TS_LIB=.agent-fuse Rscript dev/benchmarks/b3b1_endstate_probe.R

suppressMessages({
  library(TreeSearch,
          lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-fuse"),
                                  winslash = "/", mustWork = TRUE))
  library(TreeTools)
})

TARGET   <- 1261          # canonical Zanol2014 optimum (headtohead_phase0.csv tnt)
CELL_DIR <- "dev/benchmarks/basin_pull/Zanol2014"

# --- Same "-"->"?" Fitch regime the basin harness scored under ----------------
.FitchPhy <- function(p) {
  m <- PhyDatToMatrix(p, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}
loadZanol <- function() {
  e <- new.env()
  utils::data("inapplicable.phyData", package = "TreeSearch", envir = e)
  .FitchPhy(e[["inapplicable.phyData"]][["Zanol2014"]])
}

# --- Kernel data bundle (contrast / tip_data / weight / levels) ---------------
phyBundle <- function(phy) {
  at <- attributes(phy)
  list(phy = phy, contrast = at$contrast,
       tip_data = matrix(unlist(phy, use.names = FALSE),
                         nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, labels = names(phy),
       nTip = length(phy))
}

# Run the in-kernel unrooted-complete TBR (gate ON) to convergence on `tree`.
# This is the mechanism the sector gate disables; if it improves an end-state
# tree, that improvement was unreachable in the gated pipeline.
kernelUnrootedTbr <- function(tree, b) {
  edge <- Preorder(RenumberTips(tree, b$labels))[["edge"]]
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, b$contrast, b$tip_data, b$weight, b$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, unrooted = TRUE)
  out <- structure(list(edge = res$edge, Nnode = b$nTip - 1L,
                        tip.label = b$labels), class = "phylo")
  TreeLength(out, b$phy)
}

# Cheap SPR-only improving screen (O(n^2)); SPR subset of TBR, so an improving
# SPR move proves the tree is not even SPR-locally-optimal.
sprImproves <- function(tree, b) {
  base <- TreeLength(tree, b$phy)
  rr <- Preorder(RootTree(RenumberTips(tree, b$labels), b$labels[1]))
  mv <- SPRMoves(rr)
  if (!length(mv)) return(NA_real_)
  base - min(vapply(mv, TreeLength, double(1), b$phy))   # >0 => improving SPR
}

# --- Load end-state INDEP trees that missed the optimum -----------------------
phy <- loadZanol()
b <- phyBundle(phy)
cat(sprintf("Zanol2014: %d tips, %d patterns, regime '-'->'?'\n",
            b$nTip, length(b$weight)))

cells <- list.files(CELL_DIR, pattern = "^cell_.*rds$", full.names = TRUE)
parts <- lapply(cells, readRDS)
indep <- Filter(function(p) identical(p$variant, "indep"), parts)
cat(sprintf("Loaded %d cells; %d are INDEP single-replicate restarts.\n",
            length(parts), length(indep)))

# One representative end-state tree per INDEP cell (all its MPTs share `score`).
endTrees <- lapply(indep, function(p) {
  tr <- p$trees
  if (inherits(tr, "phylo")) tr else tr[[1L]]
})
storedScore <- vapply(indep, function(p) as.double(p$score), double(1))

# --- Regime-validation gate: TreeLength MUST reproduce the stored score -------
chk <- vapply(seq_along(endTrees), function(i)
  TreeLength(endTrees[[i]], phy), double(1))
drift <- abs(chk - storedScore)
cat(sprintf("Regime check: max |TreeLength - storedScore| = %.3f over %d trees\n",
            max(drift), length(drift)))
if (max(drift) > 0.5) {
  stop("REGIME MISMATCH: reloaded scoring does not reproduce stored scores; ",
       "the '-'->'?' transform or lib engine is wrong. Aborting (would be a ",
       "phantom-witness result -- see na-validation-alignment-gotcha).")
}

cat(sprintf("\nStored INDEP scores: min=%g median=%g max=%g; %d/%d reach %g.\n",
            min(storedScore), stats::median(storedScore), max(storedScore),
            sum(storedScore <= TARGET + 0.5), length(storedScore), TARGET))

missed <- which(storedScore > TARGET + 0.5)
cat(sprintf("%d INDEP trees MISS the optimum (score > %g) -- the test set.\n\n",
            length(missed), TARGET))

# --- THE TEST: unrooted-complete kernel TBR on each missed end-state tree ------
res <- lapply(missed, function(i) {
  before <- storedScore[i]
  after  <- kernelUnrootedTbr(endTrees[[i]], b)
  spr    <- sprImproves(endTrees[[i]], b)
  data.frame(cell = i, before = before, afterUnrootedTBR = after,
             tbrGain = before - after, sprGain = spr)
})
res <- do.call(rbind, res)

improved <- res[res$tbrGain > 0.5, , drop = FALSE]
cat("=== RESULT ===\n")
cat(sprintf("Missed end-state trees tested: %d\n", nrow(res)))
cat(sprintf("Improved by unrooted-complete TBR: %d (%.0f%%)\n",
            nrow(improved), 100 * nrow(improved) / nrow(res)))
if (nrow(improved)) {
  cat(sprintf("  gain: median=%g max=%g; reaches optimum (%g) on %d\n",
              stats::median(improved$tbrGain), max(improved$tbrGain), TARGET,
              sum(improved$afterUnrootedTBR <= TARGET + 0.5)))
  cat(sprintf("  of these, SPR-reachable (sprGain>0.5): %d (rest are reroot/TBR-specific)\n",
              sum(improved$sprGain > 0.5, na.rm = TRUE)))
}
cat("\nPer-tree:\n")
print(res, row.names = FALSE)

saveRDS(res, "dev/benchmarks/b3b1_endstate_result.rds")
cat("\nsaved -> dev/benchmarks/b3b1_endstate_result.rds\n")

if (nrow(improved) == 0L) {
  cat("\n=> VERDICT: 0-improving. The missed trees are already TRUE unrooted-TBR\n",
      "   optima -- the sector unrooted-completeness gate is NOT costing final\n",
      "   quality on Zanol.  B3xB1 REFUTED as a mission lever.  The residual is\n",
      "   basin/budget (matches B1 + KPI), not completeness.  Pivot to B2.\n", sep = "")
} else {
  cat("\n=> VERDICT: >0-improving. Final quality sits BELOW the unrooted-TBR local\n",
      "   optimum on some misses -- completeness IS leaving quality on the table.\n",
      "   Next: attribute (does a trailing plain-mode TBR run gate-ON already? is\n",
      "   the win SPR-reachable or reroot-specific?) before any sector_mask remap.\n", sep = "")
}
