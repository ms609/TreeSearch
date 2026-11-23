#!/usr/bin/env Rscript
# Mission A precursor -- ratchetPerturbMaxMoves 5 -> auto/scaled CORPUS ANYTIME test.
#
# WHY: the reweighting kick `ratchetPerturbMaxMoves` moved 5432's floor 1949->1945,
# but the fixed default 5 is well below the engine auto (max(20,min(200,n/8))) for
# every tree >=~40 tips, and the reach benefit was 5432-ONLY (n=1). The OPEN question
# (memory project5432-basin-structure: "wall benefit on solved matrices UNTESTED") is
# whether the bigger kick helps -- or REGRESSES -- wall-clock time-to-optimum broadly.
# The regression risk is on SMALL/EASY matrices: auto=20 is a near-random-restart on a
# 25-tip tree, which could waste re-descent and slow time-to-optimum.
#
# DESIGN (advisor-reviewed):
#  * BUILDLESS. Drive the real production `strategy="auto"` presets; override ONLY
#    ratchetPerturbMaxMoves (top-level ... arg is preserved over the preset -- the
#    fixed merge path). So this measures the REAL deployable change.
#  * Three arms in ONE roll (settles both candidate fixes named in memory):
#      fixed5  = 5   (current production default)
#      auto    = 0   (engine auto = max(20, min(200, n_tip/8)); floors at 20 for n<=160)
#      scaled  = max(5, round(n_tip/8))  (the "scale with n_tip" candidate; NO 20-floor,
#                so it stays gentle on small trees -- diverges hard from `auto` exactly
#                in the small-tree risk zone: n=25 -> scaled 5 vs auto 20)
#  * METRIC = per-replicate best-vs-wall ANYTIME trace via progressCallback (improvement
#    events only -> the step function). nThreads=1 => deterministic + serial
#    candidates_evaluated + per-replicate callback (parallel disables hits tracking).
#  * PRODUCTION STOPPING: strategy=auto, maxReplicates=96 (production), targetHits=auto.
#    So total-wall-to-stop is the real cost (a costlier kick inflates the confirm tail);
#    time-to-first-hit-of-best is the pure anytime axis. Both are read post-hoc.
#  * VALIDATION SPLIT SEQUESTERED: training matrices only (asserts split=="training").
#  * REGIME = EW Fitch, gaps->missing (the 5432 regime; `-`->`?` then MatrixToPhyDat).
#
# One SLURM array task = one (matrix, seed) cell, runs all 3 arms on the same node
# (fair same-CPU wall comparison). TASK_ID (1-based) selects the manifest row.
#
# Env: TS_LIB (install lib), NEOTRANS_DIR (matrices), CAT_CSV (catalogue),
#      OUT_DIR (output), TASK_ID / SLURM_ARRAY_TASK_ID, N_SEEDS (default 5),
#      TS_MAXREP (default 96), SMOKE (if set: 1 matrix, tiny budget -- local check).

suppressMessages({
  ts_lib <- Sys.getenv("TS_LIB", "")
  if (nzchar(ts_lib)) {
    library(TreeSearch, lib.loc = normalizePath(ts_lib, winslash = "/", mustWork = TRUE))
  } else {
    library(TreeSearch)
  }
  library(TreeTools)
})

# ---- Paths ----
find_first_dir <- function(cands) {
  for (d in cands) if (nzchar(d) && dir.exists(d)) return(normalizePath(d, winslash = "/"))
  NA_character_
}
neo_dir <- find_first_dir(c(
  Sys.getenv("NEOTRANS_DIR", ""),
  file.path(getwd(), "..", "neotrans", "inst", "matrices"),
  "C:/Users/pjjg18/GitHub/neotrans/inst/matrices",
  file.path("/nobackup", Sys.getenv("USER"), "neotrans", "inst", "matrices")
))
if (is.na(neo_dir)) stop("neotrans matrices dir not found; set NEOTRANS_DIR")

cat_csv <- Sys.getenv("CAT_CSV", "")
if (!nzchar(cat_csv) || !file.exists(cat_csv)) {
  for (p in c(file.path(getwd(), "dev/benchmarks/mbank_catalogue.csv"),
              file.path(getwd(), "mbank_catalogue.csv"))) {
    if (file.exists(p)) { cat_csv <- p; break }
  }
}
if (!file.exists(cat_csv)) stop("mbank_catalogue.csv not found; set CAT_CSV")

out_dir <- Sys.getenv("OUT_DIR", "dev/benchmarks/kick_anytime_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Fixed training sample (copied from bench_datasets.R MBANK_FIXED_SAMPLE) ----
MBANK_FIXED_SAMPLE <- c(
  # Small (20-30 taxa)
  "project532", "project2346", "project2451", "project4501",
  "project944", "project971_(1)", "project2762",
  # Medium (31-60 taxa)
  "project826", "project561", "project571", "project4146_(3)",
  "project3688", "project4049", "project423",
  # Large (61-120 taxa)
  "project4286", "project4359", "project4397", "project2084_(1)",
  "project2771", "project2184", "project3938",
  # XLarge (121+ taxa)
  "syab07201", "project4133", "project804", "project4284"
)

catalogue <- read.csv(cat_csv, stringsAsFactors = FALSE)
rownames(catalogue) <- catalogue$key

# ---- EW Fitch, gaps->missing (the 5432 regime) ----
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

# ---- Arms: value of ratchetPerturbMaxMoves to pass, per tree size ----
arm_kick <- function(arm, nTip) {
  switch(arm,
    fixed5 = 5L,
    auto   = 0L,                                   # engine: max(20, min(200, n/8))
    scaled = as.integer(max(5L, round(nTip / 8))), # no 20-floor
    stop("unknown arm: ", arm))
}
# What the engine ACTUALLY uses (for documentation/analysis).
arm_effective <- function(arm, nTip) {
  switch(arm,
    fixed5 = 5L,
    auto   = as.integer(max(20L, min(200L, round(nTip / 8)))),
    scaled = as.integer(max(5L, round(nTip / 8))))
}
ARMS <- c("fixed5", "auto", "scaled")

# ---- Anytime tracer: record best-score improvement events (step function) ----
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
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- One arm on one (matrix, seed) ----
run_arm <- function(pd, nTip, arm, seed, maxrep, cap_s) {
  kick <- arm_kick(arm, nTip)
  set.seed(seed)
  t0 <- proc.time()["elapsed"]
  tr <- make_tracer(t0)
  res <- suppressWarnings(MaximizeParsimony(
    pd, strategy = "auto",
    maxReplicates = maxrep, maxSeconds = cap_s,
    nThreads = 1L, verbosity = 0L,
    ratchetPerturbMaxMoves = kick,
    progressCallback = tr$cb))
  wall <- as.double(proc.time()["elapsed"] - t0)
  final_score <- as.double(attr(res, "score"))
  reps <- attr(res, "replicates"); reps <- if (is.null(reps)) NA_integer_ else as.integer(reps)
  cand <- attr(res, "candidates_evaluated"); cand <- if (is.null(cand)) NA_real_ else as.double(cand)
  rep_scores <- attr(res, "replicate_scores")
  trace <- if (length(tr$env$rows)) do.call(rbind, tr$env$rows) else
    data.frame(replicate = NA_integer_, elapsed_s = NA_real_,
               engine_elapsed = NA_real_, best_score = final_score)
  list(kick = kick, eff = arm_effective(arm, nTip), wall = wall,
       final_score = final_score, reps = reps, cand = cand,
       trace = trace, n_events = length(tr$env$rows),
       rep_scores_len = length(rep_scores))
}

# ---- Manifest: (matrix, seed) cells ----
N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "5"))
BASE_SEED <- 3847L
SMOKE <- nzchar(Sys.getenv("SMOKE", ""))
# maxReplicates: 300 (>= production 96, so no rep-starvation vs production) bounded
# by the tier maxSeconds cap. time-to-first-hit is read from the trace regardless of
# when the search stops, so running generously only ADDS anytime-curve, never less.
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
# maxSeconds cap per tier (safety bound; targetHits/convergence usually stops sooner).
# The kick divergence shows in the first ~minute of replicates, so caps are kept
# modest -- enough to reach the optimum and see the anytime curve, x3 arms per task
# stays well under the SLURM 1:30 limit even for the 4062-tip outlier.
cap_s <- if (SMOKE) 20 else switch(as.character(tier),
  small = 60, medium = 120, large = 240, xlarge = 480)

cat(sprintf("=== kick_anytime task %d/%d: %s (%dt, %dc, %s) seed=%d cap=%gs maxrep=%d ===\n",
            tid, nrow(manifest), key, nTip, nChar, tier, seed, cap_s, maxrep))

pd <- load_matrix(key)
stopifnot(length(pd) == nTip)

all_rows <- list()
for (arm in ARMS) {
  r <- run_arm(pd, nTip, arm, seed, maxrep, cap_s)
  cat(sprintf("  %-7s kick=%2d(eff=%3d)  final=%.0f  reps=%s  wall=%.1fs  events=%d  candM=%.2f\n",
              arm, r$kick, r$eff, r$final_score,
              ifelse(is.na(r$reps), "?", as.character(r$reps)),
              r$wall, r$n_events, r$cand / 1e6))
  # SELF-CHECK: a completed search MUST emit >=1 improvement event (rep-1 from Inf).
  if (r$n_events == 0L)
    cat(sprintf("  WARN %s: progressCallback emitted 0 improvement events -- trace facility may be broken!\n", arm))
  tr <- r$trace
  all_rows[[length(all_rows) + 1L]] <- data.frame(
    dataset = key, nTip = nTip, nChar = nChar, tier = as.character(tier),
    seed = seed, arm = arm, kick_moves = r$kick, kick_effective = r$eff,
    event = "improve", replicate = tr$replicate, elapsed_s = tr$elapsed_s,
    engine_elapsed = tr$engine_elapsed, best_score = tr$best_score,
    final_score = r$final_score, reps_done = r$reps, wall_total_s = r$wall,
    candidates = r$cand, stringsAsFactors = FALSE)
}

D <- do.call(rbind, all_rows)
of <- file.path(out_dir, sprintf("cell_%03d_%s_s%d.csv", tid, gsub("[^A-Za-z0-9]", "", key), seed))
write.csv(D, of, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", of, nrow(D)))
