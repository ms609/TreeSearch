# NA-certification gate: one cell = (dataset x seed), all arms on the same node.
#
# QUESTION.  `exact_verify_sweep` certifies that an apparently converged NA tree
# really is an unrooted-TBR optimum, and it is 97.7% of tbr_search wall when it
# runs (dev/profiling/na-exact-verify-dominates.md).  TBRParams::certify_unrooted
# lets callers whose tree is never reported skip it.  Does that pay?
#
# WHY THIS IS NOT A WALL BENCHMARK.  Under Brazeau's three-pass the indirect scan
# is only approximate, so the sweep is the ONLY exact step on the NA path -- it
# does not merely prove optimality, it finds improvers the scan cannot see.  A
# local probe on Zanol2014 (74t, default preset, 2 replicates) measured 1321
# certified vs 1326 not: five steps, for 40 s.  So this is a quality-for-speed
# trade and the primary metric is FLOOR ATTAINMENT -- the share of runs reaching
# the best score any arm found for that matrix.  Never an
# improvements-per-second ratio: one such metric once rated a 0%-reach arm a
# 1.57x win (reach-not-escape-rate).
#
# ARMS.  Three at a FIXED REPLICATE budget, which exposes the trade directly
# (quality on the left, wall on the right):
#   A_certify   today: TS_NA_NOCERTIFY unset, every caller certifies
#   B_gate      TS_NA_NOCERTIFY=1: internal sub-searches skip certification
#   C_final     B_gate + TS_NA_FINAL_CERTIFY=1: one certification of the tree
#               each replicate contributes to the pool
# ...then two more at arm A's OWN measured wall, because a wall win can be spent
# on replicates and the fixed-replicate comparison would miss that:
#   B_gate_mw / C_final_mw   maxSeconds = wall(A_certify), maxReplicates = 1000
# Matched-wall arms are re-run inside the same cell (same node, same seed), so
# the budget is genuinely comparable rather than borrowed from another node.
#
# tabuSize is a FIRST-CLASS ARM DIMENSION, not a detail.  do_reroot -- the gate
# on the certifier -- requires tabu_size == 0, and the shipped presets set
# tabuSize = 100 (default) / 200 (thorough).  At tabuSize = 100 only the sector
# sub-searches and the fuse cleanup reach the certifier; at tabuSize = 0
# (the `sprint` preset, and every raw ts_tbr_search / ts_ratchet_search call,
# which is how the 97.7% was measured) every whole-tree search does.  Running
# both is what makes the result generalise beyond one preset.
#
# `naDiag` is recorded per run.  A null result with n_evs_skipped == 0 means the
# flag never fired, which is a different finding from "it fired and bought
# nothing" -- do not let the panel conflate them.
#
# Cell index: arg[1] or $SLURM_ARRAY_TASK_ID (0-based) into
# expand.grid(dataset, seed, tabu).
# Env: TS_LIB, TS_DATASETS, TS_SEEDS, TS_REPS, TS_TABU, TS_CONCAVITY, PARTIAL_DIR.
# Local test:
#   TS_REPS=2 TS_DATASETS=Vinther2008 TS_SEEDS=1 TS_TABU=100 \
#     Rscript dev/benchmarks/bench_na_certify_cell.R 0

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                             winslash = "/"))
  library(TreeTools)
})

args  <- commandArgs(trailingOnly = TRUE)
idx   <- as.integer(if (length(args) >= 1L) args[[1]] else
                    Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
reps  <- as.integer(Sys.getenv("TS_REPS", "8"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5")),
                             "\\s+")[[1]])
tabus <- as.integer(strsplit(trimws(Sys.getenv("TS_TABU", "100 0")),
                             "\\s+")[[1]])
conc  <- as.double(Sys.getenv("TS_CONCAVITY", "Inf"))
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_na_certify")

data("inapplicable.phyData", package = "TreeSearch")
# Default = every bundled inapplicable matrix (all 30 carry "-").  The profiling
# note used four; four matrices cannot support a paired test (a sign test at
# n = 4 bottoms out at p = 0.125), and the unit of replication here is the
# MATRIX, never (matrix x seed) -- that pairing once manufactured p = 0.0007 from
# an effect that was p = 0.98 at matrix level (pair-on-matrices-not-seeds).
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
         paste(names(inapplicable.phyData), collapse = " "))), "\\s+")[[1]]

grid <- expand.grid(dataset = dsN, seed = seeds, tabu = tabus,
                    stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid)) {
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
}
row <- grid[idx + 1L, ]

d <- inapplicable.phyData[[row$dataset]]   # NATIVE NA -- do NOT recode
stopifnot(!is.null(d))

setArm <- function(noCertify, finalCertify) {
  if (noCertify) Sys.setenv(TS_NA_NOCERTIFY = "1") else Sys.unsetenv("TS_NA_NOCERTIFY")
  if (finalCertify) Sys.setenv(TS_NA_FINAL_CERTIFY = "1") else
    Sys.unsetenv("TS_NA_FINAL_CERTIFY")
}

runArm <- function(arm, noCertify, finalCertify, nReps, maxSec) {
  setArm(noCertify, finalCertify)
  set.seed(row$seed)
  # The budget -- and only the budget -- must decide when a run ends, or a faster
  # arm stops early for an unrelated reason and the arms stop being comparable.
  # stopPatience = 0 disables it.  targetHits does NOT: the C++ test is
  # `hits_to_best() >= target_hits`, so targetHits = 0 is satisfied by the first
  # replicate and truncates the run to one (measured, not assumed).  9999 is out
  # of reach at these budgets while staying small enough to avoid the huge-value
  # overflow that disabling-stopping-rules records.
  tm <- system.time(r <- suppressWarnings(MaximizeParsimony(
    d, concavity = conc, strategy = "default",
    maxReplicates = nReps, targetHits = 9999L, maxSeconds = maxSec,
    nThreads = 1L, verbosity = 0L,
    tabuSize = row$tabu, stopPatience = 0L)))
  setArm(FALSE, FALSE)
  nd <- attr(r, "naDiag")
  data.frame(
    dataset = row$dataset, seed = row$seed, tabu = row$tabu, arm = arm,
    reps_asked = nReps, maxSeconds = maxSec,
    reps_done = attr(r, "replicates"),
    score = attr(r, "score"),
    wall = as.numeric(tm[["elapsed"]]),
    nEvs = if (is.null(nd)) NA_real_ else nd$n_evs,
    nEvsSkipped = if (is.null(nd)) NA_real_ else nd$n_evs_skipped,
    nEvsImproved = if (is.null(nd)) NA_real_ else nd$n_evs_improved,
    stringsAsFactors = FALSE)
}

fixed <- rbind(
  runArm("A_certify", FALSE, FALSE, reps, 0),
  runArm("B_gate",    TRUE,  FALSE, reps, 0),
  runArm("C_final",   TRUE,  TRUE,  reps, 0)
)

# Matched wall: give the gated arms exactly the wall arm A just used, and enough
# replicates that the clock is what binds.  maxSeconds is a double, so pass A's
# measured wall unrounded -- ceiling() to whole seconds inflates a 0.7 s budget by
# 40%, which on the small matrices is larger than any effect being measured.
# maxSeconds is only polled at replicate boundaries, so the gated arms still
# overshoot by up to one replicate; the analyser prints their achieved wall so
# that overshoot is visible rather than assumed away.
mwBudget <- max(0.01, fixed$wall[fixed$arm == "A_certify"])
matched <- rbind(
  runArm("B_gate_mw",  TRUE, FALSE, 1000L, mwBudget),
  runArm("C_final_mw", TRUE, TRUE,  1000L, mwBudget)
)

out <- rbind(fixed, matched)
out$mwBudget <- mwBudget
dir.create(partdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(partdir, sprintf("nacert_%04d.csv", idx)),
          row.names = FALSE)

for (i in seq_len(nrow(out))) {
  cat(sprintf(
    "cell %d: %s seed %d tabu %d | %-11s score %.4f | %5.1f s | reps %d | evs %.0f (skip %.0f, impr %.0f)\n",
    idx, out$dataset[i], out$seed[i], out$tabu[i], out$arm[i], out$score[i],
    out$wall[i], out$reps_done[i], out$nEvs[i], out$nEvsSkipped[i],
    out$nEvsImproved[i]))
}
