#!/usr/bin/env Rscript
# Recover the TREE for the project4284 result that the ab6 confirmation A/B rests on,
# and measure the poolSuboptimal leak while we are here.
#
# WHY: reach_ab6.R recorded only attr(res, "score"). project4284 (4062 tips) is the matrix
# carrying the entire ab6 win (5/5 seeds, 353-357 vs base 359-364), and its deltas6 arm
# completed ZERO replicates -- it returned a tree found mid-replicate-0 when the deadline
# fired. A score attribute with no tree cannot be re-scored, and the standing rule here is
# to persist every new best tree the same turn. Every project4284 best_score on record in
# dev/benchmarks/ is 1040-1411 (all 30-120s budget-starved), so ~353 is the best value seen
# for this matrix anywhere -- exactly the thing that must not exist as a bare number.
#
# v2 of this script. v1 (job 18128376) died in TreeLength with "`tree` must be binary":
# MaximizeParsimony defaults to collapse = TRUE, which contracts zero-length branches, and
# TreeLength refuses a non-binary tree. Not a bug in either -- a mis-specified check.
#
# So run with collapse = FALSE, which is both re-scorable AND a free measurement. With the
# levers passed as DOTS (gate-free engine) escalatedPool is FALSE, so the return-path filter
# that normally strips poolSuboptimal trees does NOT fire -- meaning the returned set should
# contain suboptimal pool trees if poolSuboptimal = 3 really leaks. That leak is currently
# listed as unit-tested-only in reach_escalation_FINDINGS.md "Not measured"; comparing
# max(rescored) against min(rescored) measures it directly, on real data.
#
# Config is otherwise IDENTICAL to reach_ab6.R's deltas6 arm (same seeds, cap, levers,
# preprocessing), so scores should reproduce; reproducing them re-validates determinism.
# collapse is post-processing on the returned pool and does not alter the search.
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
cat(sprintf("=== recover4284 v2 task %d: seed=%d targetHits=%d cap=%gs (%d tips) ===\n",
            tid, seed, targetHits, cap_s, nTip))

set.seed(seed); t0 <- proc.time()["elapsed"]
res <- suppressWarnings(do.call(MaximizeParsimony, c(
  list(pd, strategy = "auto", maxReplicates = 300L, maxSeconds = cap_s,
       targetHits = targetHits, nThreads = 1L, verbosity = 0L,
       collapse = FALSE),                     # <- binary trees, and exposes the pool verbatim
  DELTAS6)))
wall <- as.double(proc.time()["elapsed"] - t0)
reported <- as.double(attr(res, "score"))
trees <- if (inherits(res, "phylo")) structure(list(res), class = "multiPhylo") else res
cat(sprintf("  reported=%.0f  nTrees=%d  reps=%s  wall=%.1fs\n", reported, length(trees),
            attr(res, "replicates") %||% "?", wall))

# WRITE THE TREES FIRST. v1 of this script lost a 22-minute search because it died in the
# verification step before persisting anything; the tree is the deliverable, the check is not.
allf <- file.path(out_dir, sprintf("project4284_deltas6_s%d_ALL.tre", seed))
ape::write.tree(trees, file = allf)
cat(sprintf("  wrote all %d returned tree(s) -> %s\n", length(trees), basename(allf)))

# Now verify: resolution first (TreeLength demands binary), then re-score BY LABEL.
nEdge <- vapply(trees, function(tr) nrow(tr$edge), integer(1))
binaryOK <- nEdge == 2L * nTip - 3L
cat(sprintf("  fully resolved (unrooted binary, %d edges): %d/%d   edge counts seen: %s\n",
            2L * nTip - 3L, sum(binaryOK), length(trees),
            paste(sort(unique(nEdge)), collapse = ",")))
lengths_ <- rep(NA_real_, length(trees))
for (i in seq_along(trees)) {
  lengths_[i] <- tryCatch(as.double(TreeLength(trees[[i]], pd, concavity = Inf)),
                          error = function(e) { cat(sprintf("  TreeLength failed on tree %d: %s\n",
                                                            i, conditionMessage(e))); NA_real_ })
}
ok <- !is.na(lengths_)
agree <- any(ok) && isTRUE(all.equal(min(lengths_[ok]), reported))
cat(sprintf("  re-scored %d/%d: min=%s max=%s\n", sum(ok), length(trees),
            if (any(ok)) sprintf("%.0f", min(lengths_[ok])) else "NA",
            if (any(ok)) sprintf("%.0f", max(lengths_[ok])) else "NA"))
cat(sprintf("  RE-SCORE AGREES WITH REPORTED SCORE: %s\n", agree))
if (any(ok) && !agree)
  cat("  !! MISMATCH -- reported score NOT reproduced by TreeLength on any returned tree.\n")

# The poolSuboptimal-leak measurement: any returned tree scoring worse than the best is a
# suboptimal pool tree that reached the caller, which is what the escalatedPool filter exists
# to prevent (it is inert here by design -- levers as dots, so escalatedPool == FALSE).
nAbove <- if (any(ok)) sum(lengths_[ok] > min(lengths_[ok]) + 1e-9) else NA_integer_
cat(sprintf("  SUBOPTIMAL trees in the returned set: %s of %d  -> leak %s\n",
            nAbove, sum(ok),
            if (is.na(nAbove)) "UNKNOWN" else if (nAbove > 0) "REAL (filter is load-bearing)"
            else "not observed on this cell"))

if (any(ok)) {
  best <- trees[ok][lengths_[ok] <= min(lengths_[ok]) + 1e-9]
  bf <- file.path(out_dir, sprintf("project4284_deltas6_s%d_BEST_score%.0f.tre",
                                   seed, min(lengths_[ok])))
  ape::write.tree(best, file = bf)
  cat(sprintf("  wrote %d best tree(s) -> %s\n", length(best), basename(bf)))
}
write.csv(data.frame(dataset = KEY, nTip = nTip, seed = seed, arm = "deltas6",
                     reported_score = reported,
                     rescored_min = if (any(ok)) min(lengths_[ok]) else NA_real_,
                     rescored_max = if (any(ok)) max(lengths_[ok]) else NA_real_,
                     rescore_agrees = agree, n_trees = length(trees),
                     n_rescored = sum(ok), n_suboptimal = nAbove,
                     n_binary = sum(binaryOK), all_binary = all(binaryOK),
                     reps = attr(res, "replicates") %||% NA_integer_,
                     wall_s = wall, cap_s = cap_s, collapse_arg = FALSE,
                     stringsAsFactors = FALSE),
          file.path(out_dir, sprintf("recover_%02d_s%d.csv", tid, seed)), row.names = FALSE)
cat("  CSV written.\n")
