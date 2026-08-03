#!/usr/bin/env Rscript
# T-364/T-370 battery, step 2: PROBE P -- the orientation-rate gate, and at the
# same time the behavioural discriminator that proves the three installed
# libraries really differ.
#
# WHY THIS IS A GATE.  Arm 2's pathology (T-384) fires only when the finished
# constrained Wagner tree displays the split with the tip-0-INCLUDING side as the
# rooted clade -- historically ~9.5% of starts.  Whether that ever happens
# depends on the constraint's geometry.  If the frozen constraint set happens
# never to produce it, arm 2 shows nothing, question 2 is unanswerable, and the
# whole array is wasted.  So measure the rate on the constraints we intend to
# ship, BEFORE spending cluster time.
#
# GATE: arm 2 complement-only rate materially > 0 on at least one shape;
#       arm 3 complement-only rate ~ 0; arm 1 compliance < 100%.
# All three are also the discriminator: three silently identical libraries would
# otherwise produce a clean and meaningless "no cost".
#
# Addition orders are supplied EXPLICITLY and identically to every arm, so the
# comparison is exactly paired and cannot be confounded by the arms consuming
# the RNG stream differently.
#
# Env: TS_LIB (the arm's lib), ARM (label), CONS_DIR (frozen constraints),
#      OUT_DIR, COMMON, N_ORDER (default 200), KEYS (optional subset).

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
out_dir <- Sys.getenv("OUT_DIR", "t364_probe")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# Addition orders per matrix.  Scaled DOWN with size only because the
# classification is R-level set work that grows with tip count, not because big
# matrices matter less; the rate being estimated is a property of the constraint
# geometry, which the small/medium tiers already sample densely.  xlarge is
# skipped outright (0) -- 4062 tips makes the R-side split comparison the
# dominant cost by orders of magnitude.
TIER_ORDERS <- c(small = 200L, medium = 200L, large = 50L, xlarge = 0L)
n_order_env <- Sys.getenv("N_ORDER", "")
n_order_for <- function(tier) {
  if (nzchar(n_order_env)) as.integer(n_order_env) else TIER_ORDERS[[tier]]
}

cat(sprintf("=== probe P, arm %s, TreeSearch %s ===\n",
            arm, as.character(packageVersion("TreeSearch"))))

catalogue <- load_catalogue()
files <- list.files(cons_dir, pattern = "^cons_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no frozen constraints in ", cons_dir)
keys_env <- Sys.getenv("KEYS", "")
wanted <- if (nzchar(keys_env)) trimws(strsplit(keys_env, ",")[[1]]) else NULL

rows <- list()
for (f in files) {
  ob <- readRDS(f)
  if (!is.null(wanted) && !(ob$key %in% wanted)) next
  n_order <- n_order_for(ob$tier)
  if (n_order < 1L) {
    cat(sprintf("  %-16s SKIPPED (tier %s)\n", ob$key, ob$tier))
    next
  }
  pd <- load_matrix(ob$key, catalogue)
  tips <- names(pd)
  stopifnot(identical(tips, ob$tips))
  n_tip <- length(tips)

  # Identical addition orders for every arm.
  set.seed(90210L)
  orders <- lapply(seq_len(n_order), function(i) sample.int(n_tip))

  for (nm in names(ob$cons)) {
    grp <- ob$cons[[nm]]$group
    cpd <- constraint_phydat(grp, tips)
    cc <- cx <- ok <- 0L
    walls <- numeric(n_order)
    for (i in seq_len(n_order)) {
      t0 <- proc.time()["elapsed"]
      tr <- suppressWarnings(AdditionTree(pd, constraint = cpd,
                                          sequence = orders[[i]]))
      walls[i] <- as.double(proc.time()["elapsed"] - t0)
      cl <- classify_tree(tr, grp, ob$tip0)
      cc <- cc + cl$clade_canonical
      cx <- cx + cl$clade_complement
      ok <- ok + cl$compliant
    }
    # complement-ONLY is the arm-2 signature: the tree DISPLAYS the constraint
    # but the pre-T-384 mapping cannot see it, so every regraft gets rejected.
    comp_only <- ok - cc
    cat(sprintf("  %-16s %-9s compliant=%3d/%d  canonical=%3d  complementONLY=%3d (%.1f%%)  medWall=%.3fs\n",
                ob$key, nm, ok, n_order, cc, comp_only,
                100 * comp_only / n_order, stats::median(walls)))
    rows[[length(rows) + 1L]] <- data.frame(
      arm = arm, key = ob$key, n_tip = ob$n_tip, tier = ob$tier, shape = nm,
      group_size = length(grp), has_tip0 = ob$cons[[nm]]$has_tip0,
      n_order = n_order, n_compliant = ok, n_canonical = cc,
      n_complement_only = comp_only,
      pct_violating = 100 * (n_order - ok) / n_order,
      pct_complement_only = 100 * comp_only / n_order,
      wall_median = stats::median(walls), wall_total = sum(walls),
      stringsAsFactors = FALSE)
  }
}

D <- do.call(rbind, rows)
of <- file.path(out_dir, sprintf("probe_arm%s.csv", arm))
write.csv(D, of, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", of, nrow(D)))
