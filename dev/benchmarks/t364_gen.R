#!/usr/bin/env Rscript
# T-364/T-370 battery, step 1: GENERATE AND FREEZE the constraint set.
#
# Run ONCE, against ONE arm (arm 3, the merged state).  Every arm then consumes
# the identical frozen constraints, so no arm can be advantaged by a constraint
# it helped choose.  Output: one RDS per matrix holding the enforced tip-label
# group per shape, the reference unconstrained score, and the reference tree.
#
# Two shapes per matrix, both "clade + distant rogue" (see gen_constraint),
# stratified on dataset-tip-0 membership of the enforced group:
#     rogue_no0  -- group excludes tip 0
#     rogue_in0  -- group includes tip 0
# Canonicalisation puts tip 0 outside the mask, so membership flips which side
# has to be the rooted clade: the axis that decides whether arm 2 can express
# its pathology.  Running both means the battery cannot be silently blind to it.
#
# Env: TS_LIB, NEOTRANS_DIR, CAT_CSV, OUT_DIR, COMMON (path to t364_common.R),
#      TASK_ID / SLURM_ARRAY_TASK_ID (1..25 selects the matrix), KEYS (optional
#      comma-separated override), REF_REP (reference budget, default 12).

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

out_dir <- Sys.getenv("OUT_DIR", "t364_cons")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
catalogue <- load_catalogue()

keys_env <- Sys.getenv("KEYS", "")
keys <- if (nzchar(keys_env)) trimws(strsplit(keys_env, ",")[[1]]) else MBANK_FIXED_SAMPLE
tid <- as.integer(Sys.getenv("TASK_ID", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
if (tid < 1 || tid > length(keys)) {
  stop(sprintf("TASK_ID %d out of range 1..%d", tid, length(keys)))
}
key <- keys[tid]

row <- catalogue[key, ]
n_tip <- as.integer(row$ntax)
tier <- as.character(cut(n_tip, breaks = c(0, 30, 60, 120, Inf),
                         labels = c("small", "medium", "large", "xlarge")))
ref_rep <- as.integer(Sys.getenv("REF_REP", "12"))
ref_cap <- switch(tier, small = 60, medium = 150, large = 300, xlarge = 900)

cat(sprintf("=== t364_gen %s (%dt, %dc, %s) refRep=%d cap=%gs ===\n",
            key, n_tip, as.integer(row$nchar), tier, ref_rep, ref_cap))

pd <- load_matrix(key, catalogue)
stopifnot(length(pd) == n_tip)
tips <- names(pd)
tip0 <- tips[1]

# ---- Reference UNCONSTRAINED search: supplies the tree the constraint must
# ---- conflict with, and the score the penalty is measured against.
set.seed(20260729L)
t0 <- proc.time()["elapsed"]
ref <- suppressWarnings(MaximizeParsimony(
  pd, strategy = "auto", maxReplicates = ref_rep, maxSeconds = ref_cap,
  nThreads = 1L, verbosity = 0L))
ref_wall <- as.double(proc.time()["elapsed"] - t0)
score_u <- as.double(attr(ref, "score"))
tree_u <- if (inherits(ref, "phylo")) ref else ref[[1]]
# Tip ordering must match the dataset, or "tip 0" means different things on the
# two sides of every comparison (memory na-validation-alignment-gotcha).
stopifnot(setequal(tree_u$tip.label, tips))
cat(sprintf("  unconstrained: score=%.0f  wall=%.1fs  nMPT=%d\n",
            score_u, ref_wall, if (inherits(ref, "phylo")) 1L else length(ref)))

# Deterministic tie-break helper, seeded per matrix.
set.seed(4242L + tid)
rng <- function(n) sample.int(n, 1L)

size_lo <- 3L
size_hi <- max(5L, floor(n_tip / 6))

shapes <- list(rogue_no0 = FALSE, rogue_in0 = TRUE)
cons <- list()
for (nm in names(shapes)) {
  g <- gen_constraint(tree_u, tip0, want_tip0 = shapes[[nm]],
                      size_lo = size_lo, size_hi = size_hi, rng = rng)
  if (is.null(g)) {
    cat(sprintf("  shape %-9s: NO admissible clade (size band %d-%d)\n",
                nm, size_lo, size_hi))
    next
  }
  # Sanity: the returned group must land in the requested tip-0 stratum, and
  # must be a legal split for .PrepareConstraint (1 <= |group| < nTip - 1).
  stopifnot(g$has_tip0 == shapes[[nm]],
            length(g$group) >= 1L, length(g$group) < n_tip - 1L)
  cons[[nm]] <- g
  cat(sprintf("  shape %-9s: |group|=%d (clade %d + rogue '%s')  hasTip0=%s\n",
              nm, length(g$group), g$mrca_size, g$rogue, g$has_tip0))
}
if (!length(cons)) stop("no constraint shape could be generated for ", key)

# ---- Penalty: does the constraint actually bite?  Matched budget, so this is
# ---- the operationally-relevant statement "the constrained search cannot reach
# ---- the unconstrained optimum".  It is an UPPER bound on the true optimum gap.
for (nm in names(cons)) {
  cpd <- constraint_phydat(cons[[nm]]$group, tips)
  set.seed(20260729L)
  t0 <- proc.time()["elapsed"]
  cr <- suppressWarnings(MaximizeParsimony(
    pd, constraint = cpd, strategy = "auto", maxReplicates = ref_rep,
    maxSeconds = ref_cap, nThreads = 1L, verbosity = 0L))
  cw <- as.double(proc.time()["elapsed"] - t0)
  sc <- as.double(attr(cr, "score"))
  ctree <- if (inherits(cr, "phylo")) cr else cr[[1]]
  cls <- classify_tree(ctree, cons[[nm]]$group, tip0)
  cons[[nm]]$score_constrained <- sc
  cons[[nm]]$penalty <- sc - score_u
  cons[[nm]]$ref_compliant <- cls$compliant
  cons[[nm]]$ref_wall <- cw
  cat(sprintf("  shape %-9s: constrained=%.0f  penalty=%+.0f  compliant=%s  wall=%.1fs\n",
              nm, sc, sc - score_u, cls$compliant, cw))
}

obj <- list(key = key, n_tip = n_tip, n_char = as.integer(row$nchar),
            tier = tier, tips = tips, tip0 = tip0,
            score_unconstrained = score_u, ref_rep = ref_rep, ref_cap = ref_cap,
            ref_wall = ref_wall, tree_u = tree_u, cons = cons,
            size_band = c(size_lo, size_hi))
of <- file.path(out_dir, sprintf("cons_%s.rds", gsub("[^A-Za-z0-9]", "", key)))
saveRDS(obj, of)
cat(sprintf("Wrote %s\n", of))
