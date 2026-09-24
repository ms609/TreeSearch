#!/usr/bin/env Rscript
# T-364/T-370 battery, step 3: ONE ARM of one (matrix, shape, seed) cell.
#
# The SLURM wrapper calls this three times per task -- once per arm, on the same
# node, inside one wall-clock window, in an order rotated by task id.  Three arms
# are three different INSTALLED LIBRARIES (three commits), which cannot coexist in
# one R session, hence one process per arm rather than one process per cell.
#
# WHAT IS RECORDED, and why each axis is needed:
#  * final score + total wall at a FIXED replicate budget -- the equal-work
#    comparison.  Under fixed work a blocked search may end SOONER and merely
#    score worse, so score and wall must both be reported; forcing either into
#    the headline alone would misrepresent the other.
#  * every progressCallback event (phase, elapsed, best_score, phase_score) --
#    per-replicate wall comes from consecutive replicate events, which is what
#    makes the arm-2 claim CONDITIONAL ("the cost is confined to replicates whose
#    start was complement-rooted") rather than an aggregate that averages the
#    9.5% event away.  Phase-level rows also support the T-384 signature test:
#    every phase slower on a score-identical trajectory.
#  * rooting-explicit compliance of every returned tree, so a wall difference is
#    never read without knowing whether the arms answered the same question.
#
# Env: TS_LIB, ARM, CONS_DIR, OUT_DIR, COMMON, TASK_ID/SLURM_ARRAY_TASK_ID,
#      N_SEEDS (default 3), MAXREP override, CAP_S override, MANIFEST (csv).

suppressMessages({
  ts_lib <- Sys.getenv("TS_LIB", "")
  if (nzchar(ts_lib)) {
    library(TreeSearch, lib.loc = normalizePath(ts_lib, winslash = "/", mustWork = TRUE))
  } else {
    library(TreeSearch)
  }
  library(TreeTools)
})
source(Sys.getenv("COMMON", "t364_common.R"))

arm <- Sys.getenv("ARM", "?")
cons_dir <- Sys.getenv("CONS_DIR", "t364_cons")
out_dir <- Sys.getenv("OUT_DIR", "t364_batt")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
catalogue <- load_catalogue()

# ---- Manifest: one row per (key, shape, seed_idx), built from the frozen
# ---- constraints so a matrix that yielded only one shape simply has fewer rows.
build_manifest <- function(cons_dir, n_seeds) {
  files <- sort(list.files(cons_dir, pattern = "^cons_.*\\.rds$", full.names = TRUE))
  if (!length(files)) stop("no frozen constraints in ", cons_dir)
  rows <- list()
  for (f in files) {
    ob <- readRDS(f)
    for (nm in names(ob$cons)) {
      for (s in seq_len(n_seeds)) {
        rows[[length(rows) + 1L]] <- data.frame(
          file = f, key = ob$key, shape = nm, seed_idx = s,
          n_tip = ob$n_tip, tier = ob$tier, stringsAsFactors = FALSE)
      }
    }
  }
  m <- do.call(rbind, rows)
  m[order(m$key, m$shape, m$seed_idx), ]
}

n_seeds <- as.integer(Sys.getenv("N_SEEDS", "3"))
manifest <- build_manifest(cons_dir, n_seeds)
rownames(manifest) <- NULL
mf <- Sys.getenv("MANIFEST", "")
if (nzchar(mf)) write.csv(manifest, mf, row.names = FALSE)

tid <- as.integer(Sys.getenv("TASK_ID", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
if (tid < 1) stop(sprintf("TASK_ID %d < 1", tid))
# The array is submitted at the maximum possible size (a matrix that yielded only
# one constraint shape contributes fewer rows), so overshoot is expected and must
# exit CLEAN -- a nonzero exit here would be indistinguishable from a real failure
# when the run is audited later.
if (tid > nrow(manifest)) {
  cat(sprintf("TASK_ID %d beyond manifest (%d rows) -- nothing to do\n",
              tid, nrow(manifest)))
  quit(save = "no", status = 0L)
}
mrow <- manifest[tid, ]
ob <- readRDS(mrow$file)
shape <- mrow$shape
grp <- ob$cons[[shape]]$group
BASE_SEED <- 71640L
seed <- BASE_SEED + mrow$seed_idx - 1L

# Fixed work in REPLICATES, not seconds (memory policy-in-replicates-not-seconds):
# wall is then a measurement, not a setting.  maxSeconds is a runaway guard only,
# set well above the expected cost so it does not truncate a slow arm and mask the
# very effect being measured; cells that hit it are flagged censored, never
# silently averaged in.
tier_rep <- c(small = 24L, medium = 24L, large = 24L, xlarge = 8L)
tier_cap <- c(small = 300, medium = 600, large = 1200, xlarge = 2400)
maxrep <- as.integer(Sys.getenv("MAXREP", tier_rep[[ob$tier]]))
cap_s <- as.double(Sys.getenv("CAP_S", tier_cap[[ob$tier]]))

# TIERS restricts a supplementary sweep to part of the corpus (the budget-scaling
# run at the production maxReplicates = 96 is affordable on small/medium only).
# Skipped tasks exit CLEAN, for the same auditability reason as manifest overshoot.
tiers_env <- Sys.getenv("TIERS", "")
if (nzchar(tiers_env)) {
  allowed <- trimws(strsplit(tiers_env, ",")[[1]])
  if (!(ob$tier %in% allowed)) {
    cat(sprintf("task %d: tier %s not in TIERS=%s -- skipped\n",
                tid, ob$tier, tiers_env))
    quit(save = "no", status = 0L)
  }
}

pd <- load_matrix(ob$key, catalogue)
tips <- names(pd)
stopifnot(identical(tips, ob$tips))
cpd <- constraint_phydat(grp, tips)

cat(sprintf("=== t364_cell task %d/%d arm %s | %s %s seed=%d | %dt %s rep=%d cap=%gs ===\n",
            tid, nrow(manifest), arm, ob$key, shape, seed, ob$n_tip, ob$tier,
            maxrep, cap_s))

# ---- Full event tracer: EVERY callback event, not just improvements.
`%||%` <- function(a, b) if (is.null(a)) b else a
make_tracer <- function(t0) {
  env <- new.env(parent = emptyenv())
  env$rows <- list()
  cb <- function(info) {
    env$rows[[length(env$rows) + 1L]] <- data.frame(
      phase = as.character(info$phase %||% NA_character_),
      replicate = as.integer(info$replicate %||% NA_integer_),
      wall_s = as.double(proc.time()["elapsed"] - t0),
      engine_elapsed = as.double(info$elapsed %||% NA_real_),
      best_score = as.double(info$best_score %||% NA_real_),
      phase_score = as.double(info$phase_score %||% NA_real_),
      hits_to_best = as.integer(info$hits_to_best %||% NA_integer_),
      stringsAsFactors = FALSE)
    invisible()
  }
  list(cb = cb, env = env)
}

set.seed(seed)
t0 <- proc.time()["elapsed"]
tr <- make_tracer(t0)
res <- suppressWarnings(MaximizeParsimony(
  pd, constraint = cpd, strategy = "auto",
  maxReplicates = maxrep, maxSeconds = cap_s,
  nThreads = 1L, verbosity = 0L, progressCallback = tr$cb))
wall <- as.double(proc.time()["elapsed"] - t0)

score <- as.double(attr(res, "score"))
reps <- as.integer(attr(res, "replicates") %||% NA_integer_)
cand <- as.double(attr(res, "candidates_evaluated") %||% NA_real_)
trees <- if (inherits(res, "phylo")) list(res) else res
cls <- lapply(trees, classify_tree, group = grp, tip0 = ob$tip0)
n_compliant <- sum(vapply(cls, function(x) x$compliant, logical(1)))
n_canonical <- sum(vapply(cls, function(x) x$clade_canonical, logical(1)))

# Censoring must be visible: cap_s reached means wall is a setting, not a result.
censored <- wall >= cap_s * 0.98

cat(sprintf("  score=%.0f (uncon %.0f, penalty %+.0f)  wall=%.2fs  reps=%s  nMPT=%d  compliant=%d/%d  canonical=%d  censored=%s\n",
            score, ob$score_unconstrained, score - ob$score_unconstrained, wall,
            ifelse(is.na(reps), "?", as.character(reps)), length(trees),
            n_compliant, length(trees), n_canonical, censored))

ev <- if (length(tr$env$rows)) do.call(rbind, tr$env$rows) else
  data.frame(phase = NA_character_, replicate = NA_integer_, wall_s = NA_real_,
             engine_elapsed = NA_real_, best_score = NA_real_,
             phase_score = NA_real_, hits_to_best = NA_integer_)

meta <- data.frame(
  task = tid, arm = arm, key = ob$key, n_tip = ob$n_tip, n_char = ob$n_char,
  tier = ob$tier, shape = shape, has_tip0 = ob$cons[[shape]]$has_tip0,
  group_size = length(grp), seed = seed, maxrep = maxrep, cap_s = cap_s,
  score = score, score_unconstrained = ob$score_unconstrained,
  penalty = score - ob$score_unconstrained,
  gen_penalty = ob$cons[[shape]]$penalty,
  wall_s = wall, reps_done = reps, candidates = cand, n_mpt = length(trees),
  n_compliant = n_compliant, n_canonical = n_canonical, censored = censored,
  n_events = nrow(ev), pkg_version = as.character(packageVersion("TreeSearch")),
  stringsAsFactors = FALSE)

tag <- sprintf("t%03d_arm%s", tid, arm)

# Per-PHASE cumulative wall (ms), straight from the engine.  Needed to test the
# recorded T-384 diagnostic signature -- "every phase slower on a score-identical
# trajectory" -- on its own terms.  The progressCallback only emits `replicate` and
# `done`, so without this attribute the best available test would be "the search
# ran longer", which is a much weaker claim than the recorded one and would not
# distinguish "each phase is slower" from "there are more replicates".
tm <- attr(res, "timings")
if (!is.null(tm) && length(tm)) {
  tmd <- data.frame(task = tid, arm = arm, key = ob$key, shape = shape,
                    seed = seed, phase = names(tm),
                    ms = as.double(tm), reps_done = reps, score = score,
                    stringsAsFactors = FALSE)
  write.csv(tmd, file.path(out_dir, sprintf("timings_%s.csv", tag)), row.names = FALSE)
}

write.csv(meta, file.path(out_dir, sprintf("meta_%s.csv", tag)), row.names = FALSE)
ev2 <- cbind(task = tid, arm = arm, key = ob$key, shape = shape, seed = seed, ev)
write.csv(ev2, file.path(out_dir, sprintf("events_%s.csv", tag)), row.names = FALSE)
cat(sprintf("Wrote meta_%s.csv + events_%s.csv (%d events)\n", tag, tag, nrow(ev)))
