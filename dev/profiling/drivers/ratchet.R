# Ratchet inner-loop profiling driver — Area #2
# bare: 3.7 s on 2026-05-18 (thorough × 3 reps, nThreads=1)
# Dataset: Zhu2013 (75 tips, 4 states, 43% missing)
# Strategy: "default" preset (ratchetCycles=12), maxReplicates=1, nThreads=1
#
# What this exercises: ratchet_search() in ts_ratchet.cpp, which loops
#   [save_perturb → perturb → tbr_search → restore_perturb → accept/reject]
# for n_cycles iterations around an initial TBR pass.

# Default lib.loc ".vtune-lib" preserves the local VTune convention; CI overrides
# via TREESEARCH_VTUNE_LIB (same install path, kept uniform with the other drivers).
LIBDIR <- Sys.getenv("TREESEARCH_VTUNE_LIB", unset = ".vtune-lib")
library(TreeSearch, lib.loc = LIBDIR)

set.seed(5813)

# TS_DATASET/TS_REPS default to Zhu2013/3 => byte-identical to the local VTune
# driver; CI overrides to a small dataset so callgrind (~40x) completes under cap.
ds_name <- Sys.getenv("TS_DATASET", unset = "Zhu2013")
n_reps  <- as.integer(Sys.getenv("TS_REPS", unset = "3"))
dataset <- inapplicable.phyData[[ds_name]]

# Suppress replicate-count adequacy warning (1 rep is intentional)
t0 <- proc.time()
result <- suppressWarnings(
  MaximizeParsimony(
    dataset,
    maxReplicates = n_reps,
    targetHits    = 1L,
    nThreads      = 1L,
    strategy      = "thorough",
    verbosity     = 0L
  )
)
elapsed <- round((proc.time() - t0)["elapsed"], 1)
cat("Elapsed:", elapsed, "s | Score:", attr(result, "score"), "\n")
