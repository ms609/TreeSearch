# T-373 repro: with `total_words == 0` every search kernel bails, but under
# HSJ/XFORM the objective is still topology-dependent -- so the start tree is
# returned unsearched, with no warning.
#
# Red-team area 10, 2026-07-28 (opus verifier). Ported into the repo from the
# session scratchpad because findings.md cited it by a temp path that gets
# swept. Standalone; NOT wired into the test suite.
#
#   Rscript dev/red-team/heavy-tests/hsj-totalwords-zero-noop.R [path/to/library]
#
# WHY THIS SHAPE, and why the obvious comparison is invalid: the tempting test
# is "search returns 8, but a random sample of trees contains one scoring 6".
# That is CONFOUNDED -- `RandomTree(ds, root = TRUE)` roots every draw on tip 1
# while the engine reports its score at its own internal rooting, and HSJ/XFORM
# scores are rooting-dependent (T-374). So a score gap can be pure rooting.
#
# This design avoids rooting entirely. It pins `startEdge`, turns search effort
# up to maximal, and asks a binary question: is the returned edge matrix
# BIT-IDENTICAL to the start? The `oneW` arm is the control -- one single
# non-zero Fitch weight, the SAME hierarchy DP -- and it must move the tree.
#
# Result at 1a94403b / 5cffb18d: zeroW returns the start tree bit-identical
# (10 -> 10) while oneW moves it (13 -> 8). 2250 of the 10395 8-tip topologies
# score 9 in that state, so the returned 10 is genuinely sub-optimal.
#
# CAVEAT: this calls `TreeSearch:::ts_driven_search` with an explicit argument
# list, so it is coupled to that internal signature and will need updating if it
# changes. If it errors on argument mismatch, that is script rot, not a fix.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) .libPaths(c(args[[1]], .libPaths()))
suppressMessages({library(TreeSearch); library(TreeTools)})

tips <- paste0("t", 1:8)
mat <- rbind(
  t1 = c("0", "0", "1", "1", "0", "1", "1"),
  t2 = c("0", "1", "1", "1", "1", "1", "0"),
  t3 = c("1", "0", "0", "0", "-", "0", "-"),
  t4 = c("1", "1", "0", "1", "1", "1", "1"),
  t5 = c("0", "1", "0", "0", "-", "1", "0"),
  t6 = c("1", "0", "1", "1", "0", "0", "-"),
  t7 = c("1", "1", "1", "0", "-", "1", "1"),
  t8 = c("0", "0", "0", "1", "1", "0", "-")
)
rownames(mat) <- tips
ds <- phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                       ambiguity = "?")
hierarchy <- CharacterHierarchy("4" = 5, "6" = 7)

at       <- attributes(ds)
tipData  <- matrix(unlist(ds, use.names = FALSE), nrow = length(ds), byrow = TRUE)
blocks   <- TreeSearch:::.HierarchyToBlocks(hierarchy)
tipLabs  <- TreeSearch:::.BuildTipLabels(ds)
absState <- TreeSearch:::.HSJAbsentState(ds)
hsjCfg   <- list(hierarchyBlocks = blocks, hsjTipLabels = tipLabs,
                 hsjAlpha = 1.0, hsjAbsentState = absState)

Score <- function(edge, w) {
  TreeSearch:::ts_hsj_score(edge, at$contrast, tipData, as.integer(w),
                            at$levels, blocks, 1.0, tipLabs, absState)
}

zeroW <- integer(7)                 # every Fitch weight 0 -> total_words == 0
oneW  <- c(1L, rep(0L, 6))          # total_words > 0, SAME hierarchy DP

# A deliberately poor start tree, so there is plenty of room to improve.
start <- Renumber(RenumberTips(
  ape::read.tree(text = "(t1,(t2,(t3,(t4,(t5,(t6,(t7,t8)))))));"), tips))

cat(sprintf("start tree HSJ score -- zeroW: %s   oneW: %s\n\n",
            Score(start$edge, zeroW), Score(start$edge, oneW)))

Run <- function(w, seed = 5) {
  # Maximal effort: if anything CAN rearrange the tree, it will.
  ctl <- SearchControl(tbrMaxHits = 200L, ratchetCycles = 50L,
                       driftCycles = 20L, xssRounds = 3L, rssRounds = 3L,
                       cssRounds = 3L, fuseInterval = 0L,
                       poolMaxSize = 1L, poolSuboptimal = 0)
  rt <- list(maxReplicates = 1L, targetHits = 1000L, maxSeconds = 0,
             verbosity = 0L, nThreads = 1L,
             startEdge = list(start$edge), progressCallback = NULL)
  sg <- list(min_steps = integer(0), concavity = Inf, xpiwe = FALSE,
             xpiwe_r = 0.5, xpiwe_max_f = 5, obs_count = integer(0),
             infoAmounts = NULL)
  set.seed(seed)
  TreeSearch:::ts_driven_search(at$contrast, tipData, as.integer(w), at$levels,
                               ctl, rt, sg, NULL, hsjCfg, NULL)
}

noop <- c(zeroW = NA, oneW = NA)
for (nm in c("zeroW", "oneW")) {
  w <- get(nm)
  out <- Run(w)
  same <- identical(as.integer(out$trees[[1]]), as.integer(start$edge))
  noop[[nm]] <- same
  cat(sprintf("%-6s total_words>0? %-5s | best=%-6s start=%-6s | returned == start: %s\n",
              nm, nm != "zeroW", out$best_score, Score(start$edge, w),
              if (same) "YES (no rearrangement ever applied)" else "no (search worked)"))
}

cat("\n---- verdict ----\n")
if (isTRUE(noop[["zeroW"]]) && isFALSE(noop[["oneW"]])) {
  cat("T-373 REPRODUCES: zeroW is a silent no-op while the control searches.\n")
  quit(status = 1L)
}
if (isTRUE(noop[["oneW"]])) {
  cat("INCONCLUSIVE: the control also failed to move the tree, so this run\n")
  cat("proves nothing about the zeroW arm. Investigate the harness first.\n")
  quit(status = 2L)
}
cat("zeroW searched (T-373 would be fixed).\n")
