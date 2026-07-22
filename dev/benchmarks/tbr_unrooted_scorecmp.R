# tbr_unrooted_scorecmp.R -- validate the DIRECT unrooted root-edge path against
# the PHYSICAL-REROOT reference (TS_PHYS_REROOT=1), by comparing the kernel's own
# converged score (result.best_score) on identical starts.  Physical reroot is
# complete by construction (tries all rootings, exact per-scorer scoring), so
# direct == phys on every tree => the direct path is equally complete.  Uses the
# kernel-native score, so NO TreeLength scoring-match is required (avoids the
# apples-to-oranges trap for IW/NA).
#
# Usage: Rscript dev/benchmarks/tbr_unrooted_scorecmp.R [nTrees] [nTip] [concavity]
#   concavity < 0  => equal weights (EW);  finite >0 => implied weights (IW)
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})
args   <- commandArgs(trailingOnly = TRUE)
nTrees <- if (length(args) >= 1) as.integer(args[[1]]) else 60L
nTip   <- if (length(args) >= 2) as.integer(args[[2]]) else 12L
conc   <- if (length(args) >= 3) as.numeric(args[[3]]) else -1
nChar  <- 60L; nState <- 3L

randomData <- function(seed) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  m <- matrix(sample(0:(nState - 1L), nTip * nChar, replace = TRUE),
              nrow = nTip, dimnames = list(tips, NULL))
  phy <- phangorn::phyDat(m, type = "USER", levels = as.character(0:(nState - 1L)))
  at <- attributes(phy)
  list(phy = phy, contrast = at$contrast,
       tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, nTip = length(phy), labels = names(phy))
}

runK <- function(tree, d, seed, phys) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  if (phys) Sys.setenv(TS_PHYS_REROOT = "1") else Sys.unsetenv("TS_PHYS_REROOT")
  set.seed(seed)
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L,
    concavity = conc, unrooted = TRUE)
  Sys.unsetenv("TS_PHYS_REROOT")
  res$score
}

cat(sprintf("=== direct vs physical-reroot score-cmp: %d trees, %d tips, concavity=%s (%s) ===\n",
            nTrees, nTip, conc, if (conc < 0) "EW" else "IW"))
mism <- 0L; worseDirect <- 0L
for (i in seq_len(nTrees)) {
  d <- randomData(1000L + i)
  set.seed(7000L + i); start <- RandomTree(d$phy, root = TRUE)
  sD <- runK(start, d, i, FALSE)
  sP <- runK(start, d, i, TRUE)
  if (abs(sD - sP) > 1e-6) {
    mism <- mism + 1L
    if (sD > sP + 1e-6) worseDirect <- worseDirect + 1L
    cat(sprintf("  tree %d: direct=%.4f  phys=%.4f  (%s)\n", i, sD, sP,
                if (sD > sP) "DIRECT WORSE (incomplete)" else "direct better"))
  }
}
cat(sprintf("\nMISMATCHES: %d / %d  (direct strictly worse: %d)\n", mism, nTrees, worseDirect))
if (mism == 0L) cat("=> direct path reaches the SAME optimum as physical reroot on all trees.\n")
