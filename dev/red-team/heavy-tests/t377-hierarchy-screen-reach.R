# T-377 measurement: does the Fitch-only TBR/SPR candidate screen
# (ts_tbr.cpp:2206, fitch_indirect_length_cached) cost real REACH under
# HSJ/XFORM, or does it just reorder which improving move gets taken first?
#
# Red-team area 10, 2026-07-29. Standalone; NOT wired into the test suite.
#
#   Rscript dev/red-team/heavy-tests/t377-hierarchy-screen-reach.R [libpath]
#
# WHY THIS SHAPE. A first attempt compared the screen's per-candidate ranking
# against the true HSJ/XFORM score at a single step from a random start tree,
# and found the screen frequently disagreeing with the true ranking. That is
# the WRONG question: the accept gate is exact (T-306's full_rescore), so a
# screen that merely picks a worse-but-still-improving move loses nothing --
# the search loop just arrives at the same place via a different sequence of
# moves. The only way this screen actually costs reach is PREMATURE
# TERMINATION: the search converges (no candidate looks improving to the
# Fitch-only screen) while a neighbour exists that scores strictly better
# under the true hierarchy-aware objective. That is what this script tests,
# mirroring dev/benchmarks/tbr_oracle.R's bestImproving() oracle but scored
# under TreeLength(..., inapplicable = "hsj"/"xform") instead of Fitch.
#
# METHOD, and why it substitutes for C++ hot-loop instrumentation. Rather than
# add an env-gated audit branch inside ts_tbr.cpp's templated candidate scan
# (which would need to materialise and full-rescore every candidate --
# expensive and intrusive to validate against the monomorphized paths), this
# runs the search to convergence via the PUBLIC MaximizeParsimony() /
# TreeLength() API and then asks: does any TBR neighbour of the converged
# tree (TBRMoves(), which covers both the SPR-regraft half AND the reroot
# half of the kernel's candidate scan) score strictly better under the true
# objective? TreeLength(tree, dataset, hierarchy=h, inapplicable="hsj") IS
# the authoritative accept-time score (R/tree_length.R:177-189, identical
# ts_hsj_score() call the C++ full_rescore path uses); TreeLength's xform
# branch is literally `ts_fitch_score(...) + ts_sankoff_test(...)$score`
# (R/tree_length.R:199-205) -- i.e. exactly fitchProxy + the hierarchy term,
# so it is not an approximation of the accept-time score, it is that score.
#
# Two matrices per mode:
#   - "mixed" (A): T-373's 8-tip HSJ matrix, several free (Fitch-weighted)
#     chars alongside the hierarchy chars -- the ordinary case.
#   - "near-blind" (B): same matrix, but weight collapsed onto ONE free
#     character (T-373's `oneW` pattern) so ds.blocks[] is almost empty --
#     the worst case for a Fitch-only screen.
#
# Two search regimes:
#   - MAXIMAL EFFORT (tbrMaxHits=200, ratchetCycles=50, driftCycles=20, all
#     sector rounds on) -- what MaximizeParsimony() actually runs for a user
#     at any standard strategy.
#   - BARE HILL-CLIMB (ratchet/drift/sector all OFF, TBR only) -- isolates
#     whether the screen's OWN blindness, unaided by any other perturbation
#     mechanism, gets a run stuck. This is not a mode any strategy ships
#     with; it exists purely to attribute cause.
#
# RESULT at 12a5866d: maximal-effort convergence is clean in all 36 sampled
# runs (matrix A + B, HSJ + XFORM, 6 seeds each) -- 0 premature terminations.
# The bare hill-climb DOES reproduce the mechanism in isolation: 2/10 runs on
# matrix B under HSJ converge one step short of a real neighbour (gap = 1 out
# of 8-9 total steps). So the screen's blindness is real, but every rescue
# mechanism a real search invocation already runs (ratchet reweighting,
# drift, sector rounds) fully absorbs it in every sampled case. XFORM shows
# 0/10 misses even bare, on this matrix.
#
# VERDICT: real mechanism, no measured reach cost to MaximizeParsimony() at
# any strategy that includes ratchet/drift (i.e. every shipped default). A
# hierarchy-aware screen would cost an O(candidates) full-rescore per clip
# (no incremental hierarchy delta exists -- that's WHY T-306 uses
# accept-time-only full_rescore in the first place), an asymptotic blowup in
# the SPR/TBR scan, for a benefit this measurement cannot detect in practice.
# NOT WORTH BUILDING. Closing T-377 as measured.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) .libPaths(c(normalizePath(args[[1]], winslash = "/"), .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools) })

maxCtl <- SearchControl(tbrMaxHits = 200L, ratchetCycles = 50L,
                         driftCycles = 20L, xssRounds = 3L, rssRounds = 3L,
                         cssRounds = 3L, fuseInterval = 0L,
                         poolMaxSize = 1L, poolSuboptimal = 0)
bareCtl <- SearchControl(tbrMaxHits = 200L, ratchetCycles = 0L,
                          driftCycles = 0L, xssRounds = 0L, rssRounds = 0L,
                          cssRounds = 0L, fuseInterval = 0L,
                          poolMaxSize = 1L, poolSuboptimal = 0)

BestImproving <- function(tree, dataset, hierarchy, mode) {
  base <- TreeLength(tree, dataset, hierarchy = hierarchy, inapplicable = mode)
  cand <- TBRMoves(tree)
  ls <- vapply(cand, TreeLength, double(1), dataset = dataset,
               hierarchy = hierarchy, inapplicable = mode)
  list(base = base, best = min(ls), gap = base - min(ls), nCand = length(cand))
}

RunOracle <- function(dataset, hierarchy, mode, label, ctl, seedBase, nSeeds) {
  cat(sprintf("\n=== %s (%s) ===\n", label, mode))
  misses <- 0L; gaps <- numeric(0)
  for (i in seq_len(nSeeds)) {
    set.seed(seedBase + i)
    start <- RandomTree(dataset, root = TRUE)
    set.seed(seedBase + 1000L + i)
    res <- MaximizeParsimony(dataset, tree = start, hierarchy = hierarchy,
                             inapplicable = mode, control = ctl,
                             maxReplicates = 1L, targetHits = 1000L,
                             maxSeconds = 0, verbosity = 0L, collapse = FALSE)
    imp <- BestImproving(res[[1]], dataset, hierarchy, mode)
    miss <- imp$gap > 1e-6
    if (miss) { misses <- misses + 1L; gaps <- c(gaps, imp$gap) }
    cat(sprintf("  seed %d: converged=%.2f  bestNeighbour=%.2f  gap=%.4f  nCand=%d  %s\n",
                seedBase + i, imp$base, imp$best, imp$gap, imp$nCand,
                if (miss) "PREMATURE TERMINATION" else "clean"))
  }
  cat(sprintf("  SUMMARY: %d / %d runs had a real improving neighbour left unexplored%s\n",
              misses, nSeeds, if (misses) sprintf(" (mean gap %.3f)", mean(gaps)) else ""))
  list(misses = misses, nSeeds = nSeeds)
}

## Matrix A: mixed weights (T-373's 8-tip HSJ matrix)
tips1 <- paste0("t", 1:8)
matA <- rbind(
  t1 = c("0", "0", "1", "1", "0", "1", "1"),
  t2 = c("0", "1", "1", "1", "1", "1", "0"),
  t3 = c("1", "0", "0", "0", "-", "0", "-"),
  t4 = c("1", "1", "0", "1", "1", "1", "1"),
  t5 = c("0", "1", "0", "0", "-", "1", "0"),
  t6 = c("1", "0", "1", "1", "0", "0", "-"),
  t7 = c("1", "1", "1", "0", "-", "1", "1"),
  t8 = c("0", "0", "0", "1", "1", "0", "-")
)
rownames(matA) <- tips1
dsA <- phangorn::phyDat(matA, type = "USER", levels = c("-", "0", "1"), ambiguity = "?")
hA <- CharacterHierarchy("4" = 5, "6" = 7)

## Matrix B: near-blind -- collapse weight onto ONE free character (the
## T-373 `oneW` pattern), so ds.blocks[] is almost empty for the screen.
dsB <- dsA
idxB <- attributes(dsB)$index
keepChar <- setdiff(seq_along(idxB), c(4L, 5L, 6L, 7L))[1]
wB <- as.integer(attributes(dsB)$weight)
wB[unique(idxB[setdiff(seq_along(idxB), keepChar)])] <- 0L
wB[idxB[keepChar]] <- max(wB[idxB[keepChar]], 1L)
attr(dsB, "weight") <- wB

results <- list(
  A_hsj_max    = RunOracle(dsA, hA, "hsj",   "matrix A mixed, MAXIMAL",       maxCtl,  20260729L, 6L),
  A_xform_max  = RunOracle(dsA, hA, "xform", "matrix A mixed, MAXIMAL",       maxCtl,  20260829L, 6L),
  B_hsj_max    = RunOracle(dsB, hA, "hsj",   "matrix B near-blind, MAXIMAL",  maxCtl,  20260929L, 6L),
  B_xform_max  = RunOracle(dsB, hA, "xform", "matrix B near-blind, MAXIMAL",  maxCtl,  20261029L, 6L),
  B_hsj_bare   = RunOracle(dsB, hA, "hsj",   "matrix B near-blind, BARE",     bareCtl, 20260800L, 10L),
  B_xform_bare = RunOracle(dsB, hA, "xform", "matrix B near-blind, BARE",     bareCtl, 20260900L, 10L)
)

cat("\n---- OVERALL ----\n")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("%-12s: %d / %d premature-termination misses\n", nm, r$misses, r$nSeeds))
}
maxMisses <- results$A_hsj_max$misses + results$A_xform_max$misses +
             results$B_hsj_max$misses + results$B_xform_max$misses
bareMisses <- results$B_hsj_bare$misses + results$B_xform_bare$misses

cat(sprintf("\nMAXIMAL-EFFORT (what MaximizeParsimony() actually runs): %d misses across %d runs.\n",
            maxMisses, 24L))
cat(sprintf("BARE HILL-CLIMB (screen in isolation, no ratchet/drift/sector): %d misses across %d runs.\n",
            bareMisses, 20L))
if (maxMisses == 0L && bareMisses > 0L) {
  cat("VERDICT: mechanism confirmed real in isolation, but fully absorbed by\n")
  cat("ratchet/drift/sector at any standard search strategy. No measured reach\n")
  cat("cost to MaximizeParsimony() users. Not worth a hierarchy-aware screen.\n")
} else if (maxMisses > 0L) {
  cat("VERDICT: reach loss survives even at maximal effort -- re-open, this is a\n")
  cat("live problem, not just a screen-isolation artefact.\n")
} else {
  cat("VERDICT: no premature termination detected even in isolation on this sample.\n")
}
