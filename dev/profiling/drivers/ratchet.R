# Ratchet inner-loop profiling driver — Area #2
# bare: 3.7 s on 2026-05-18 (thorough × 3 reps, nThreads=1)
# Dataset: Zhu2013 (75 tips, 4 states, 43% missing)
# Strategy: "default" preset (ratchetCycles=12), maxReplicates=1, nThreads=1
#
# What this exercises: ratchet_search() in ts_ratchet.cpp, which loops
#   [save_perturb → perturb → tbr_search → restore_perturb → accept/reject]
# for n_cycles iterations around an initial TBR pass.

library(TreeSearch, lib.loc = ".vtune-lib")

set.seed(5813)

dataset <- inapplicable.phyData[["Zhu2013"]]

# Suppress replicate-count adequacy warning (1 rep is intentional)
t0 <- proc.time()
result <- suppressWarnings(
  MaximizeParsimony(
    dataset,
    maxReplicates = 3L,
    targetHits    = 1L,
    nThreads      = 1L,
    strategy      = "thorough",
    verbosity     = 0L
  )
)
elapsed <- round((proc.time() - t0)["elapsed"], 1)
cat("Elapsed:", elapsed, "s | Score:", attr(result, "score"), "\n")
