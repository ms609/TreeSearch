#!/usr/bin/env Rscript
# Confirmation A/B of the SHIPPED escalation form (lead 5 / "belt-and-braces").
#
# The ship-gate A/B measured SEVEN levers, UNGATED, via strategy="auto" over all 25
# MBANK_FIXED_SAMPLE matrices. What actually ships is SIX levers (ratchetCycles dropped to
# .IwRatchetDepth) GATED to strategy thorough/large. That configuration was never run as
# such; the ship argument is a tier decomposition of the earlier run. This closes that gap.
#
# Scope: only the matrices where auto resolves to thorough (61-120 tips) or large (>=121),
# i.e. the 11 large+xlarge members of MBANK_FIXED_SAMPLE -- on small/medium the gate now
# does nothing, so those cells would be identical in both arms and carry no information.
# 11 matrices x 5 seeds = 55 cells. Everything else matches the original harness: both arms
# at the SAME raised targetHits (the trigger threshold), differing ONLY in the six levers;
# gate-free engine with the levers as dots, so the result does not depend on the gate
# implementation; union-best target; training split asserted.
suppressMessages({
  ts_lib <- Sys.getenv("TS_LIB", "")
  if (nzchar(ts_lib)) {
    library(TreeSearch, lib.loc = normalizePath(ts_lib, winslash = "/", mustWork = TRUE))
  } else library(TreeSearch)
  library(TreeTools)
})
neo_dir <- Sys.getenv("NEOTRANS_DIR"); cat_csv <- Sys.getenv("CAT_CSV")
out_dir <- Sys.getenv("OUT_DIR", "."); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# large (61-120) + xlarge (>=121) members of MBANK_FIXED_SAMPLE
KEYS <- c("project4286", "project4359", "project4397", "project2084_(1)",
          "project2771", "project2184", "project3938",
          "syab07201", "project4133", "project804", "project4284")

# The SHIPPED six (ratchetCycles deliberately absent -- owned by .IwRatchetDepth).
DELTAS6 <- list(ratchetPerturbMaxMoves = 0L, driftCycles = 25L,
                postRatchetSectorial = TRUE, stallEscalateFactor = 1.5,
                intraFuse = TRUE, poolSuboptimal = 3)
reach_threshold <- function(nTip) 2L * max(10L, as.integer(nTip / 5))

catalogue <- read.csv(cat_csv, stringsAsFactors = FALSE)
rownames(catalogue) <- catalogue$key
to_fitch <- function(pd) { m <- PhyDatToMatrix(pd, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
load_matrix <- function(key) {
  row <- catalogue[key, ]
  if (!identical(row$split, "training"))
    stop(sprintf("key %s is split='%s' -- validation is SEQUESTERED", key, row$split))
  to_fitch(suppressWarnings(TreeTools::ReadAsPhyDat(file.path(neo_dir, row$filename))))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
make_tracer <- function(t0) {
  env <- new.env(parent = emptyenv()); env$prev <- Inf; env$rows <- list()
  cb <- function(info) {
    if (!identical(info$phase, "replicate")) return(invisible())
    bs <- info$best_score
    if (is.null(bs) || length(bs) != 1L || !is.finite(bs)) return(invisible())
    if (bs < env$prev - 1e-9) { env$prev <- bs
      env$rows[[length(env$rows) + 1L]] <- data.frame(
        replicate = as.integer(info$replicate %||% NA_integer_),
        elapsed_s = as.double(proc.time()["elapsed"] - t0),
        best_score = as.double(bs), stringsAsFactors = FALSE) }
    invisible()
  }
  list(cb = cb, env = env)
}
run_arm <- function(pd, nTip, arm, seed, maxrep, cap_s) {
  set.seed(seed); t0 <- proc.time()["elapsed"]; tr <- make_tracer(t0)
  args <- list(pd, strategy = "auto", maxReplicates = maxrep, maxSeconds = cap_s,
               targetHits = reach_threshold(nTip), nThreads = 1L, verbosity = 0L,
               progressCallback = tr$cb)
  if (identical(arm, "deltas6")) args <- c(args, DELTAS6)
  res <- suppressWarnings(do.call(MaximizeParsimony, args))
  wall <- as.double(proc.time()["elapsed"] - t0)
  trace <- if (length(tr$env$rows)) do.call(rbind, tr$env$rows) else
    data.frame(replicate = NA_integer_, elapsed_s = NA_real_,
               best_score = as.double(attr(res, "score")))
  list(wall = wall, final_score = as.double(attr(res, "score")),
       reps = attr(res, "replicates") %||% NA_integer_,
       cand = attr(res, "candidates_evaluated") %||% NA_real_,
       trace = trace, n_events = length(tr$env$rows))
}
N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "5")); BASE_SEED <- 7731L
maxrep <- as.integer(Sys.getenv("TS_MAXREP", "300"))
manifest <- expand.grid(key = KEYS, seed_idx = seq_len(N_SEEDS), stringsAsFactors = FALSE)
manifest <- manifest[order(manifest$key, manifest$seed_idx), ]; rownames(manifest) <- NULL
tid <- as.integer(Sys.getenv("TASK_ID", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
if (tid < 1 || tid > nrow(manifest)) stop(sprintf("TASK_ID %d out of 1..%d", tid, nrow(manifest)))
key <- manifest$key[tid]; seed <- BASE_SEED + manifest$seed_idx[tid] - 1L
row <- catalogue[key, ]; nTip <- as.integer(row$ntax)
tier <- if (nTip >= 121L) "xlarge" else "large"
cap_s <- if (identical(tier, "xlarge")) 1440 else 720
cat(sprintf("=== ab6 task %d/%d: %s (%dt, %s) seed=%d targetHits=%d cap=%gs ===\n",
            tid, nrow(manifest), key, nTip, tier, seed, reach_threshold(nTip), cap_s))
pd <- load_matrix(key); stopifnot(length(pd) == nTip)
rows <- list()
for (arm in c("base", "deltas6")) {
  r <- run_arm(pd, nTip, arm, seed, maxrep, cap_s)
  cat(sprintf("  %-8s final=%.0f reps=%s wall=%.1fs events=%d\n",
              arm, r$final_score, r$reps, r$wall, r$n_events))
  if (r$n_events == 0L) cat("  WARN: 0 improvement events -- trace may be broken\n")
  tr <- r$trace
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = key, nTip = nTip, tier = tier, seed = seed, arm = arm,
    target_hits = reach_threshold(nTip), event = "improve", replicate = tr$replicate,
    elapsed_s = tr$elapsed_s, best_score = tr$best_score, final_score = r$final_score,
    reps_done = r$reps, wall_total_s = r$wall, candidates = r$cand, cap_s = cap_s,
    stringsAsFactors = FALSE)
}
D <- do.call(rbind, rows)
of <- file.path(out_dir, sprintf("cell_%03d_%s_s%d.csv", tid, gsub("[^A-Za-z0-9]", "", key), seed))
write.csv(D, of, row.names = FALSE); cat(sprintf("Wrote %s (%d rows)\n", of, nrow(D)))
