# NA-certification under the SHIPPED STOPPING RULES.  One cell = (dataset x seed).
#
# WHY A SECOND PANEL.  The first panel (bench_na_certify_cell.R) deliberately
# neutralised targetHits and stopPatience so that the budget -- and only the
# budget -- ended a run, and it measured two budgets: fixed 8 replicates (the gate
# LOSES reach, p = 0.035) and fixed wall (the gate WINS, +0.167, p = 0.012).
# Production is NEITHER.  `MaximizeParsimony` defaults to maxSeconds = 0,
# maxReplicates = 96 and targetHits = max(10, ntax/5), so a real run stops when
# the best score has been hit targetHits times.  Nothing measured so far says
# where that regime falls, and it is the only regime that decides a preset
# default.
#
# THE MECHANISM AT ISSUE.  `targetHits` counts hits against the CURRENT best, not
# the true optimum.  A uniformly weaker (gated) search therefore has two opposed
# effects: it accumulates hits more slowly, so it runs MORE replicates and gets
# more chances (good, and the first panel says more replicates is where the gate
# wins) -- but if it reliably plateaus one step above the optimum it accumulates
# targetHits hits on THAT score and stops early, confidently wrong (bad).  Which
# effect dominates is exactly what this panel measures.
#
# ARMS.
#   A_default      today, every default untouched
#   B_gate         TS_NA_NOCERTIFY=1, every default untouched
#   C_final        B_gate + TS_NA_FINAL_CERTIFY=1 (never needed >33 replicates in
#                  panel 1, so the shipped maxReplicates = 96 cap never binds)
#   B_gate_reps    B_gate + maxReplicates = 480.  Panel 1's WINNING arm spent
#                  101-226 replicates on 23 of 30 matrices, so the 96 cap would
#                  have bound on all 23 -- gating without raising the cap is NOT
#                  the arm that won, and this is the pairing to test.
#   A_hits3        today, targetHits = 3x default
#   B_gate_hits3   B_gate,  targetHits = 3x default
# The last two answer "does anything need adjusting as targetHits rises?"
# directly: raising it is already read as a user signal that the matrix is hard
# (it scales the IW ratchet via .IwRatchetDepth), and it should push both arms
# toward the many-replicate regime where the gate won.
#
# Records reps / hitsToBest as well as score and wall, because "where did it
# stop, and why" is the whole question -- a null on score with a large reps
# difference is a finding, not a wash.
#
# Cell index: arg[1] or $SLURM_ARRAY_TASK_ID (0-based) into expand.grid(dataset, seed).
# Env: TS_LIB, TS_DATASETS, TS_SEEDS, TS_CONCAVITY, PARTIAL_DIR.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                             winslash = "/"))
  library(TreeTools)
})

args  <- commandArgs(trailingOnly = TRUE)
idx   <- as.integer(if (length(args) >= 1L) args[[1]] else
                    Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5")),
                             "\\s+")[[1]])
conc  <- as.double(Sys.getenv("TS_CONCAVITY", "Inf"))
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_na_certify_stop")

data("inapplicable.phyData", package = "TreeSearch")
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
         paste(names(inapplicable.phyData), collapse = " "))), "\\s+")[[1]]

grid <- expand.grid(dataset = dsN, seed = seeds, stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid)) {
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
}
row <- grid[idx + 1L, ]
d <- inapplicable.phyData[[row$dataset]]
stopifnot(!is.null(d))

# Mirrors MaximizeParsimony's own default: max(10L, as.integer(NTip / 5)).
defaultHits <- max(10L, as.integer(length(d) / 5))

runArm <- function(arm, noCertify, finalCertify, extra = list()) {
  if (noCertify) Sys.setenv(TS_NA_NOCERTIFY = "1") else Sys.unsetenv("TS_NA_NOCERTIFY")
  if (finalCertify) Sys.setenv(TS_NA_FINAL_CERTIFY = "1") else
    Sys.unsetenv("TS_NA_FINAL_CERTIFY")
  set.seed(row$seed)
  # NOTE: no targetHits / stopPatience / maxSeconds overrides.  Every stopping
  # rule is left exactly as shipped -- that is the point of this panel.
  cl <- c(list(d, concavity = conc, strategy = "default",
               nThreads = 1L, verbosity = 0L), extra)
  tm <- system.time(r <- suppressWarnings(do.call(MaximizeParsimony, cl)))
  Sys.unsetenv("TS_NA_NOCERTIFY"); Sys.unsetenv("TS_NA_FINAL_CERTIFY")
  nd <- attr(r, "naDiag")
  data.frame(
    dataset = row$dataset, ntax = length(d), seed = row$seed, arm = arm,
    defaultHits = defaultHits,
    targetHits = if (is.null(extra$targetHits)) defaultHits else extra$targetHits,
    maxReplicates = if (is.null(extra$maxReplicates)) 96L else extra$maxReplicates,
    score = attr(r, "score"),
    reps = attr(r, "replicates"),
    hitsToBest = attr(r, "hits_to_best"),
    timedOut = isTRUE(attr(r, "timed_out")),
    perturbStop = isTRUE(attr(r, "perturb_stop")),
    wall = as.numeric(tm[["elapsed"]]),
    nEvs = if (is.null(nd)) NA_real_ else nd$n_evs,
    nEvsSkipped = if (is.null(nd)) NA_real_ else nd$n_evs_skipped,
    nEvsImproved = if (is.null(nd)) NA_real_ else nd$n_evs_improved,
    stringsAsFactors = FALSE)
}

out <- rbind(
  runArm("A_default",    FALSE, FALSE),
  runArm("B_gate",       TRUE,  FALSE),
  runArm("C_final",      TRUE,  TRUE),
  runArm("B_gate_reps",  TRUE,  FALSE, list(maxReplicates = 480L)),
  runArm("A_hits3",      FALSE, FALSE, list(targetHits = 3L * defaultHits)),
  runArm("B_gate_hits3", TRUE,  FALSE, list(targetHits = 3L * defaultHits))
)

dir.create(partdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(partdir, sprintf("nastop_%04d.csv", idx)),
          row.names = FALSE)
for (i in seq_len(nrow(out))) {
  cat(sprintf(
    "cell %d: %s(%dt) seed %d | %-13s score %.4f | reps %3d hits %3d (target %d) | %6.1f s | evs %.0f skip %.0f\n",
    idx, out$dataset[i], out$ntax[i], out$seed[i], out$arm[i], out$score[i],
    out$reps[i], out$hitsToBest[i], out$targetHits[i], out$wall[i],
    out$nEvs[i], out$nEvsSkipped[i]))
}
