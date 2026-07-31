#!/usr/bin/env Rscript
# v1 REACH-ESCALATION GATE -- general-pool anytime A/B (the SHIP/NO-SHIP gate).
#
# WHAT IS BEING DECIDED: `MaximizeParsimony()` gains a gate that, when the user raises
# `targetHits` to >= 2 * max(10, nTip/5), applies 7 deeper per-replicate perturbation
# deltas (ratchetCycles 40, kick 0/auto-deep, driftCycles 25, postRatchetSectorial,
# stallEscalateFactor 1.5, intraFuse, poolSuboptimal 3). On the hard-reach tail
# (project5432) these measurably improve the score distribution. This A/B asks the
# OTHER question: when the gate fires on ORDINARY data, does it help or HURT?
#
# DESIGN (advisor-reviewed):
#  * BUILDLESS + GATE-FREE ENGINE. Runs on plain cpp-search (curlib) and passes the
#    deltas as top-level dots -- byte-equivalent to what the gate applies, so no gated
#    build is needed and the result is not contingent on my patch being correct.
#  * THE CONFOUND THIS AVOIDS: comparing default-targetHits vs raised-targetHits would
#    confound the deltas with "raised targetHits searches longer anyway". So BOTH arms
#    use the SAME raised targetHits (exactly the gate threshold) and differ ONLY in the
#    7 deltas. This isolates the escalation.
#      base   = strategy auto, targetHits = 2*max(10, nTip/5), stock preset
#      deltas = identical + the 7 escalation deltas   <- what the gate does when it fires
#  * strategy = "auto" (NOT pinned): the gate layers on whatever preset auto picks, so
#    auto is the faithful test -- and sprint/default (driftCycles 0-2) is exactly where
#    driftCycles=25 is the largest relative change, i.e. where a regression would show.
#  * METRIC = per-replicate improvement trace via progressCallback (the step function).
#    Report REACH + tt_hit (wall) + **rep2hit (replicates)**. rep2hit is mandatory: the
#    2026-07-14 kick study's apparent 7-12% wall2hit "regression" was overturned by
#    rep2hit = 1.000 -- same replicate to optimum, just costlier per replicate. Reporting
#    wall alone manufactures a false regression verdict.
#  * TARGET = union-best final score across arms within each (matrix, seed) cell -- the
#    established mbank convention (there is NO canonical best-known table for mbank).
#  * BUDGET: tier caps 3x the kick_anytime caps. driftCycles 2->25 is ~10x the drift
#    work; reusing the old caps would starve the deltas arm and guarantee a spurious
#    regression -- exactly the trap that time-truncated 5432 arm B (it stopped at
#    maxSeconds*(1-enumTimeFraction), never spending its replicate budget).
#  * VALIDATION SPLIT SEQUESTERED: asserts split == "training" (project_id %% 5 == 0 is
#    validation and is a one-way door).
#  * REGIME = EW Fitch, gaps->missing.
#
# PRE-REGISTERED DECISION RULE (fixed BEFORE any result was seen):
#   The gate is OPT-IN and fires only when the user has said "time is no issue" by
#   raising targetHits. So a wall cost is ACCEPTED by construction; what is NOT
#   acceptable is the escalation finding WORSE trees.
#     SHIP    : reach(deltas) >= reach(base) overall AND no size tier shows a reach
#               regression. Wall cost reported honestly alongside.
#     NO-SHIP : reach(deltas) < reach(base) on the general pool (the deltas actively
#               hurt when they fire) -> restrict the gate or drop v1.
#   Evidence format = paired per-cell better/worse/tie counts + median paired ratios
#   (the repo convention), NOT p-values.
#
# One SLURM array task = one (matrix, seed) cell, both arms on the same node (fair
# same-CPU wall comparison). TASK_ID (1-based) selects the manifest row.
#
# Env: TS_LIB, NEOTRANS_DIR, CAT_CSV, OUT_DIR, TASK_ID / SLURM_ARRAY_TASK_ID,
#      N_SEEDS (default 5), TS_MAXREP (default 300), SMOKE (1 matrix, tiny budget).

suppressMessages({
  ts_lib <- Sys.getenv("TS_LIB", "")
  if (nzchar(ts_lib)) {
    library(TreeSearch, lib.loc = normalizePath(ts_lib, winslash = "/", mustWork = TRUE))
  } else {
    library(TreeSearch)
  }
  library(TreeTools)
})

neo_dir <- Sys.getenv("NEOTRANS_DIR", "")
cat_csv <- Sys.getenv("CAT_CSV", "")
out_dir <- Sys.getenv("OUT_DIR", ".")
if (!nzchar(neo_dir)) stop("NEOTRANS_DIR unset")
if (!nzchar(cat_csv)) stop("CAT_CSV unset")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Fixed 25-matrix training sample (bench_datasets.R MBANK_FIXED_SAMPLE).
# "Do not modify: results are only comparable when the same sample is used."
MBANK_FIXED_SAMPLE <- c(
  "project532", "project2346", "project2451", "project4501",
  "project944", "project971_(1)", "project2762",
  "project826", "project561", "project571", "project4146_(3)",
  "project3688", "project4049", "project423",
  "project4286", "project4359", "project4397", "project2084_(1)",
  "project2771", "project2184", "project3938",
  "syab07201", "project4133", "project804", "project4284"
)

catalogue <- read.csv(cat_csv, stringsAsFactors = FALSE)
rownames(catalogue) <- catalogue$key

to_fitch <- function(pd) {
  m <- PhyDatToMatrix(pd, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}

load_matrix <- function(key) {
  if (!key %in% catalogue$key) stop("key not in catalogue: ", key)
  row <- catalogue[key, ]
  # SEQUESTER: refuse validation-split matrices.
  if (!identical(row$split, "training"))
    stop(sprintf("key %s is split='%s' -- validation is SEQUESTERED", key, row$split))
  f <- file.path(neo_dir, row$filename)
  if (!file.exists(f)) stop("matrix file not found: ", f)
  pd <- suppressWarnings(TreeTools::ReadAsPhyDat(f))
  to_fitch(pd)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# The 7 escalation deltas -- must stay byte-identical to .ReachEscalationDeltas().
ESCALATION_DELTAS <- list(
  ratchetCycles = 40L, ratchetPerturbMaxMoves = 0L, driftCycles = 25L,
  postRatchetSectorial = TRUE, stallEscalateFactor = 1.5,
  intraFuse = TRUE, poolSuboptimal = 3
)
# The gate's trigger, mirroring .ReachEscalationThreshold(nTip).
reach_threshold <- function(nTip) 2L * max(10L, as.integer(nTip / 5))
ARMS <- c("base", "deltas")

# ---- Anytime tracer: improvement events only (the step function) ----
make_tracer <- function(t0) {
  env <- new.env(parent = emptyenv())
  env$prev <- Inf
  env$rows <- list()
  cb <- function(info) {
    if (!identical(info$phase, "replicate")) return(invisible())
    bs <- info$best_score
    if (is.null(bs) || length(bs) != 1L || !is.finite(bs)) return(invisible())
    if (bs < env$prev - 1e-9) {
      env$prev <- bs
      env$rows[[length(env$rows) + 1L]] <- data.frame(
        replicate = as.integer(info$replicate %||% NA_integer_),
        elapsed_s = as.double(proc.time()["elapsed"] - t0),
        engine_elapsed = if (!is.null(info$elapsed)) as.double(info$elapsed) else NA_real_,
        best_score = as.double(bs), stringsAsFactors = FALSE)
    }
    invisible()
  }
  list(cb = cb, env = env)
}

# The engine reserves a fraction of maxSeconds for MPT enumeration and stops the main search
# at maxSeconds * (1 - ENUM_TIME_FRACTION). Passed EXPLICITLY (at the engine default, so
# behaviour is unchanged) and recorded per row, so the analyzer can compute the real deadline
# instead of guessing at the nominal cap -- guessing is what let 59/110 deadline-bound ab6
# cells be scored as "converged".
ENUM_TIME_FRACTION <- 0.1

run_arm <- function(pd, nTip, arm, seed, maxrep, cap_s) {
  set.seed(seed)
  t0 <- proc.time()["elapsed"]
  tr <- make_tracer(t0)
  args <- list(pd, strategy = "auto",
               maxReplicates = maxrep, maxSeconds = cap_s,
               enumTimeFraction = ENUM_TIME_FRACTION,
               targetHits = reach_threshold(nTip),   # SAME in both arms
               nThreads = 1L, verbosity = 0L,
               progressCallback = tr$cb)
  if (identical(arm, "deltas")) args <- c(args, ESCALATION_DELTAS)
  res <- suppressWarnings(do.call(MaximizeParsimony, args))
  wall <- as.double(proc.time()["elapsed"] - t0)
  final_score <- as.double(attr(res, "score"))
  reps <- attr(res, "replicates"); reps <- if (is.null(reps)) NA_integer_ else as.integer(reps)
  cand <- attr(res, "candidates_evaluated"); cand <- if (is.null(cand)) NA_real_ else as.double(cand)
  trace <- if (length(tr$env$rows)) do.call(rbind, tr$env$rows) else
    data.frame(replicate = NA_integer_, elapsed_s = NA_real_,
               engine_elapsed = NA_real_, best_score = final_score)
  list(wall = wall, final_score = final_score, reps = reps, cand = cand,
       trace = trace, n_events = length(tr$env$rows))
}

N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "5"))
BASE_SEED <- 5821L
SMOKE <- nzchar(Sys.getenv("SMOKE", ""))
maxrep <- as.integer(Sys.getenv("TS_MAXREP", if (SMOKE) "3" else "300"))

keys <- MBANK_FIXED_SAMPLE
if (SMOKE) { keys <- keys[1]; N_SEEDS <- 1L }

manifest <- expand.grid(key = keys, seed_idx = seq_len(N_SEEDS),
                        stringsAsFactors = FALSE)
manifest <- manifest[order(manifest$key, manifest$seed_idx), ]
rownames(manifest) <- NULL

tid <- as.integer(Sys.getenv("TASK_ID", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
if (SMOKE) tid <- 1L
if (tid < 1 || tid > nrow(manifest))
  stop(sprintf("TASK_ID %d out of range 1..%d", tid, nrow(manifest)))

key <- manifest$key[tid]
seed <- BASE_SEED + manifest$seed_idx[tid] - 1L
row <- catalogue[key, ]
nTip <- as.integer(row$ntax)
nChar <- as.integer(row$nchar)
tier <- cut(nTip, breaks = c(0, 30, 60, 120, Inf),
            labels = c("small", "medium", "large", "xlarge"))
# 3x the kick_anytime caps: driftCycles 2->25 is ~10x the drift work, so the old caps
# would starve the deltas arm and manufacture a spurious regression.
cap_s <- if (SMOKE) 20 else switch(as.character(tier),
  small = 180, medium = 360, large = 720, xlarge = 1440)

cat(sprintf("=== reach_ab task %d/%d: %s (%dt, %dc, %s) seed=%d targetHits=%d cap=%gs maxrep=%d ===\n",
            tid, nrow(manifest), key, nTip, nChar, tier, seed,
            reach_threshold(nTip), cap_s, maxrep))

pd <- load_matrix(key)
stopifnot(length(pd) == nTip)

all_rows <- list()
for (arm in ARMS) {
  r <- run_arm(pd, nTip, arm, seed, maxrep, cap_s)
  cat(sprintf("  %-6s final=%.0f  reps=%s  wall=%.1fs  events=%d  candM=%.2f\n",
              arm, r$final_score,
              ifelse(is.na(r$reps), "?", as.character(r$reps)),
              r$wall, r$n_events, r$cand / 1e6))
  # SELF-CHECK: a completed search MUST emit >=1 improvement event (rep 1 from Inf).
  if (r$n_events == 0L)
    cat(sprintf("  WARN %s: progressCallback emitted 0 improvement events -- trace facility may be broken!\n", arm))
  tr <- r$trace
  all_rows[[length(all_rows) + 1L]] <- data.frame(
    dataset = key, nTip = nTip, nChar = nChar, tier = as.character(tier),
    seed = seed, arm = arm, target_hits = reach_threshold(nTip),
    event = "improve", replicate = tr$replicate, elapsed_s = tr$elapsed_s,
    engine_elapsed = tr$engine_elapsed, best_score = tr$best_score,
    final_score = r$final_score, reps_done = r$reps, wall_total_s = r$wall,
    candidates = r$cand, cap_s = cap_s,
    enum_time_fraction = ENUM_TIME_FRACTION, stringsAsFactors = FALSE)
}

D <- do.call(rbind, all_rows)
of <- file.path(out_dir, sprintf("cell_%03d_%s_s%d.csv", tid, gsub("[^A-Za-z0-9]", "", key), seed))
write.csv(D, of, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", of, nrow(D)))
