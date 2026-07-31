#!/usr/bin/env Rscript
# Recover the TREE for the project4284 result that the ab6 confirmation A/B rests on.
#
# WHY: reach_ab6.R recorded only attr(res, "score"). project4284 (4062 tips) is the matrix
# carrying the entire ab6 win (5/5 seeds, 353-357 vs base 359-364), and its deltas6 arm
# completed ZERO replicates -- it returned a tree found mid-replicate-0 when the deadline
# fired. A score attribute with no tree cannot be re-scored, and the standing rule here is
# to persist every new best tree the same turn. Every project4284 best_score on record in
# dev/benchmarks/ is 1040-1411 (all 30-120s budget-starved), so ~353 is the best value seen
# for this matrix anywhere -- exactly the thing that must not exist as a bare number.
#
# Config is IDENTICAL to reach_ab6.R's deltas6 arm (same seeds, same cap, same levers, same
# preprocessing) so the scores should reproduce; reproducing them also re-validates
# determinism. The only change is that the tree is written out and independently re-scored
# BY LABEL (never by raw edge+tip_data index -- RenumberTips permutes).
suppressMessages({
  ts_lib <- Sys.getenv("TS_LIB", "")
  if (nzchar(ts_lib)) {
    library(TreeSearch, lib.loc = normalizePath(ts_lib, winslash = "/", mustWork = TRUE))
  } else library(TreeSearch)
  library(TreeTools)
})
neo_dir <- Sys.getenv("NEOTRANS_DIR"); cat_csv <- Sys.getenv("CAT_CSV")
out_dir <- Sys.getenv("OUT_DIR", "."); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

KEY <- "project4284"
`%||%` <- function(a, b) if (is.null(a)) b else a
DELTAS6 <- list(ratchetPerturbMaxMoves = 0L, driftCycles = 25L,
                postRatchetSectorial = TRUE, stallEscalateFactor = 1.5,
                intraFuse = TRUE, poolSuboptimal = 3)

catalogue <- read.csv(cat_csv, stringsAsFactors = FALSE)
rownames(catalogue) <- catalogue$key
row <- catalogue[KEY, ]
if (!identical(row$split, "training"))
  stop(sprintf("key %s is split='%s' -- validation is SEQUESTERED", KEY, row$split))
to_fitch <- function(pd) {
  m <- PhyDatToMatrix(pd, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m)
}
pd <- to_fitch(suppressWarnings(TreeTools::ReadAsPhyDat(file.path(neo_dir, row$filename))))
nTip <- as.integer(row$ntax); stopifnot(length(pd) == nTip)

tid <- as.integer(Sys.getenv("TASK_ID", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
seed <- 7731L + tid - 1L                      # ab6 used BASE_SEED 7731 + seed_idx - 1
cap_s <- 1440                                 # ab6 xlarge cap, unchanged
targetHits <- 2L * max(10L, as.integer(nTip / 5))
cat(sprintf("=== recover4284 task %d: seed=%d targetHits=%d cap=%gs (%d tips) ===\n",
            tid, seed, targetHits, cap_s, nTip))

set.seed(seed); t0 <- proc.time()["elapsed"]
res <- suppressWarnings(do.call(MaximizeParsimony, c(
  list(pd, strategy = "auto", maxReplicates = 300L, maxSeconds = cap_s,
       targetHits = targetHits, nThreads = 1L, verbosity = 0L), DELTAS6)))
wall <- as.double(proc.time()["elapsed"] - t0)
reported <- as.double(attr(res, "score"))

# INDEPENDENT re-score, by label, of every returned tree. TreeLength() re-derives the score
# from the tree + dataset, so agreement with attr(res,"score") is a real check that the
# returned object is a valid tree scoring what the harness claimed -- the point of the run.
trees <- if (inherits(res, "phylo")) structure(list(res), class = "multiPhylo") else res
lengths_ <- vapply(trees, function(tr) as.double(TreeLength(tr, pd, concavity = Inf)),
                   double(1))
cat(sprintf("  reported=%.0f  nTrees=%d  re-scored: min=%.0f max=%.0f  wall=%.1fs\n",
            reported, length(trees), min(lengths_), max(lengths_), wall))
agree <- isTRUE(all.equal(min(lengths_), reported))
cat(sprintf("  RE-SCORE AGREES WITH REPORTED SCORE: %s\n", agree))
if (!agree)
  cat("  !! MISMATCH -- the reported score is NOT reproduced by TreeLength on the tree.\n")
binaryOK <- vapply(trees, function(tr) length(tr$edge[, 1]) == 2L * nTip - 3L, logical(1))
cat(sprintf("  fully resolved (unrooted binary): %d/%d\n", sum(binaryOK), length(trees)))

best <- trees[lengths_ <= min(lengths_) + 1e-9]
tf <- file.path(out_dir, sprintf("project4284_deltas6_s%d_score%.0f.tre", seed, min(lengths_)))
ape::write.tree(best, file = tf)
cat(sprintf("  WROTE %d tree(s) -> %s\n", length(best), tf))
write.csv(data.frame(dataset = KEY, nTip = nTip, seed = seed, arm = "deltas6",
                     reported_score = reported, rescored_min = min(lengths_),
                     rescore_agrees = agree, n_trees = length(trees),
                     n_best = length(best), all_binary = all(binaryOK),
                     reps = attr(res, "replicates") %||% NA_integer_,
                     wall_s = wall, cap_s = cap_s, tree_file = basename(tf),
                     stringsAsFactors = FALSE),
          file.path(out_dir, sprintf("recover_%02d_s%d.csv", tid, seed)), row.names = FALSE)
