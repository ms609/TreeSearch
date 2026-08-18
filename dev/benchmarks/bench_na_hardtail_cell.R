# Does more effort recover a better score on TOUGH inapplicable matrices?
# One cell = (dataset x seed); every arm runs on the same node.
#
# WHY A THIRD PANEL.  Panels 1 and 2 answered the wrong question about effort.
# Panel 2's `A_hits3` arm raised `targetHits` x3 and improved nothing -- but it
# could not have: on the three hardest matrices the 96-replicate cap bound BOTH
# arms on every seed, so it performed byte-identical work to the control, and
# 94% of the extra replicates it bought corpus-wide landed on matrices that were
# already solved.  `targetHits` cannot act once `maxReplicates` binds, and on
# hard data it always binds first.
#
# So the effort hypothesis -- a low-effort run misses the optimum, a high-effort
# run takes longer and recovers a better score -- has never been tested.
#
# ARMS ARE `effort` ITSELF, not a hand-rolled budget.  A first draft of this
# panel raised only `maxReplicates` and was inert on Aria2015: `targetHits`
# stopped every run at 39 replicates, far below even the 96 cap, so all three
# budgets returned identical scores.  That is the exact MIRROR of panel 2's
# mistake -- where `targetHits` binds, raising the cap alone does nothing, just
# as raising `targetHits` does nothing where the cap binds.  Only moving BOTH
# escalates every dataset, which is precisely what the ladder does at rungs 5+.
# Testing the shipped argument is therefore both more honest and more
# informative than testing a knob no user turns.
#
#   effort 0 / +1 / +2, crossed with certification on/off.
#
# The certified column tests the effort hypothesis directly.  The gated column
# asks the follow-on: can the ~17x wall that gating frees be spent on more search
# to beat certification?  Panel 2 says no at 5x the replicates (Zanol2014:
# certified@96 = 1311 beats gated@480 = 1313); this is the fair rematch.
#
# NB panel 2 pinned `strategy = "default"` for every matrix, so on the 65-119-tip
# matrices it did NOT measure what a user gets -- `auto` selects `thorough`
# there.  `effort = 0` here is the genuine shipped default.
#
# TWO OTHER DESIGN CHANGES, both forced by the same mistake:
#
# 1. HARD MATRICES ONLY.  24 of the 30 bundled inapplicable matrices sit at 1.0
#    in every arm under the shipped budget; including them can only dilute.  The
#    six here are the only ones where any arm ever failed to reach the best score
#    found.
# 2. SCORE, NOT BINARY FLOOR ATTAINMENT.  On a matrix where no arm reaches a
#    known floor, attainment is 0 everywhere and blind to an arm that got
#    closer.  Panel 2 hid exactly that.  Score is the metric; attainment relative
#    to the best any arm found is reported alongside, not instead.
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
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS",
           paste(1:10, collapse = " "))), "\\s+")[[1]])
conc  <- as.double(Sys.getenv("TS_CONCAVITY", "Inf"))
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_na_hardtail")

# The six matrices where ANY arm of panel 2 fell short of the best score found.
# Zanol2014 and Aguado2009 are the only two where the shipped default itself
# falls short; the other four are marginal and kept as a sanity band.
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
  "Zanol2014 Aguado2009 Zhu2013 Wortley2006 Geisler2001 Aria2015")), "\\s+")[[1]]

data("inapplicable.phyData", package = "TreeSearch")
grid <- expand.grid(dataset = dsN, seed = seeds, stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid)) {
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
}
row <- grid[idx + 1L, ]
d <- inapplicable.phyData[[row$dataset]]
stopifnot(!is.null(d))

runArm <- function(arm, noCertify, effort) {
  if (noCertify) Sys.setenv(TS_NA_NOCERTIFY = "1") else Sys.unsetenv("TS_NA_NOCERTIFY")
  set.seed(row$seed)
  # Nothing pinned but `effort`: the point is what a user actually gets when they
  # turn the one dial up.  `stopPatience = 0` only so that a dry spell cannot end
  # a run before its budget, which would confound the arms.
  tm <- system.time(r <- suppressWarnings(MaximizeParsimony(
    d, concavity = conc, effort = effort, stopPatience = 0L,
    nThreads = 1L, verbosity = 0L)))
  Sys.unsetenv("TS_NA_NOCERTIFY")
  nd <- attr(r, "naDiag")
  data.frame(
    dataset = row$dataset, ntax = length(d), seed = row$seed, arm = arm,
    certify = !noCertify, effort = effort,
    score = attr(r, "score"),
    reps = attr(r, "replicates"),
    hitsToBest = attr(r, "hits_to_best"),
    nTopologies = attr(r, "n_topologies"),
    wall = as.numeric(tm[["elapsed"]]),
    nEvs = if (is.null(nd)) NA_real_ else nd$n_evs,
    nEvsSkipped = if (is.null(nd)) NA_real_ else nd$n_evs_skipped,
    nEvsImproved = if (is.null(nd)) NA_real_ else nd$n_evs_improved,
    stringsAsFactors = FALSE)
}

out <- rbind(
  runArm("certify_e0", FALSE, 0L),
  runArm("certify_e1", FALSE, 1L),
  runArm("certify_e2", FALSE, 2L),
  runArm("gate_e0",    TRUE,  0L),
  runArm("gate_e1",    TRUE,  1L),
  runArm("gate_e2",    TRUE,  2L)
)

dir.create(partdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(partdir, sprintf("nahard_%04d.csv", idx)),
          row.names = FALSE)
for (i in seq_len(nrow(out))) {
  cat(sprintf(
    "cell %d: %s(%dt) seed %2d | %-11s score %.4f | reps %4d hits %3d | %8.1f s\n",
    idx, out$dataset[i], out$ntax[i], out$seed[i], out$arm[i], out$score[i],
    out$reps[i], out$hitsToBest[i], out$wall[i]))
}
