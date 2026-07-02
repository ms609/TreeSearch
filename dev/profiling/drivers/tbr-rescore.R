# TBR full-rescore profiling driver — Area #4
# bare: 3.9 s on 2026-05-19 (Zhu2013 75t, 12 ratchet reps × 12 cycles)
# Dataset: Zhu2013 (75 tips, inapplicable characters)
# Strategy: ts_ratchet_search() called directly (no MaximizeParsimony overhead)
#            to match the ratchet context from which the 28 % estimate came.
#            Each ratchet cycle: perturb → tbr_search → restore → accept/reject.
#            full_rescore at ts_tbr.cpp:1138 fires on every accepted TBR move.
#
# Why ts_ratchet_search rather than ts_tbr_search directly?
#   - TBR from a near-optimal tree converges instantly (0 accepts → 0
#     full_rescore-at-acceptance events; not representative).
#   - Ratchet perturbs the tree, driving TBR from sub-optimal states that
#     generate many accepts and thus many full_rescore calls.
#   - Perturbation overhead < 2 % (confirmed in prior round); ratchet time is
#     effectively all TBR time.
#
# VTune attribution note:
#   full_rescore() is a 2-line static inline under -O2.
#   Attribution falls to reset_states() / score_tree() at source lines:
#     ts_tbr.cpp:1138  (rescore after acceptance — the T-300 target)
#     ts_tbr.cpp:563   (initial full_rescore at tbr_search entry)
#     ts_tbr.cpp:1283  (trailing full_rescore at exit)
#   Use: vtune -report hotspots -group-by source-line -filter module=TreeSearch.dll

# Timestamped lib from build 2026-05-19 06:10:49
LIBDIR <- Sys.getenv("TREESEARCH_VTUNE_LIB",
                     unset = "dev/profiling/.vtune-lib-20260519061049")
library(TreeSearch, lib.loc = LIBDIR)

set.seed(5813)

# TS_DATASET/TS_REPS default to Zhu2013/12 => byte-identical to the local VTune
# driver; CI overrides to a small dataset so callgrind (~40x) completes under cap.
ds_name <- Sys.getenv("TS_DATASET", unset = "Zhu2013")
dataset <- inapplicable.phyData[[ds_name]]
at       <- attributes(dataset)
contrast <- at$contrast
tip_data <- matrix(unlist(dataset, use.names = FALSE),
                   nrow = length(dataset), byrow = TRUE)
weight   <- TreeSearch:::.ScaleWeight(at$weight)
levels   <- at$levels

# Starting tree: random unrooted tree (cold start, like a new replicate)
starting_edge <- ape::rtree(length(dataset), tip.label = names(dataset),
                             rooted = FALSE)
starting_edge <- ape::root(starting_edge, 1L, resolve.root = TRUE)[["edge"]]
stopifnot(starting_edge[1L, 1L] > length(dataset))  # internal node first

N_REPS <- as.integer(Sys.getenv("TS_REPS", unset = "12"))
t0 <- proc.time()
for (rep in seq_len(N_REPS)) {
  set.seed(rep)
  result <- TreeSearch:::ts_ratchet_search(
    edge        = starting_edge,
    contrast    = contrast,
    tip_data    = tip_data,
    weight      = weight,
    levels      = levels,
    nCycles     = 12L,
    perturbProb = 0.04,
    maxHits     = 1L
  )
  # Keep each rep independent: restart from the same cold tree
  # (varying seed ensures different topology visits and thus more accept events)
}
elapsed <- round((proc.time() - t0)["elapsed"], 1)
cat("Elapsed:", elapsed, "s |", N_REPS, "ratchet reps (nCycles=12, Zhu2013 75t) | score:", result$score, "\n")
