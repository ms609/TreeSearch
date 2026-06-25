# ============================================================================
# Hamilton sweep: tune "thorough" to capture BOTH Zhu2013 MPT islands reliably,
# and settle wagnerStarts 3-vs-5 (the only difference between thorough and the
# opt-in "intensive" preset) so the two can be merged.
# ============================================================================
#
# See dev/benchmarks/bench_thorough_two_islands.R and data/README.md for the
# two-island pathology. The disconnected 16-tree island is reached only when a
# main-loop replicate's random start lands in its (narrow) basin; "thorough"
# should make that near-certain.
#
# Two parts, one CSV per array shard (concatenate afterwards):
#   PART A (completeness): vendored two-island fixture x VARIANTS x BUDGETS x
#     SEEDS. Metric = did the run recover BOTH islands; how many of the 16; and
#     how many replicates completed (to test the throughput hypothesis directly).
#   PART B (score/throughput regression + merge): standard benchmark datasets x
#     VARIANTS x BUDGETS x SEEDS. Metric = score-over-target, replicates, wall.
#     Confirms the completeness lever does not regress wall-clock-to-optimum and
#     re-examines the wagnerStarts trade-off behind the thorough/intensive split.
#
# VARIANTS are dots-overrides on strategy="thorough" (dots override the preset;
# the proven pattern from hamilton_thorough_rasstarts.R). Validation-set note:
# Part B uses the classic inapplicable.phyData benchmarks (the set the original
# wagnerStarts decision was made on); the mbank VALIDATION split stays sequestered.
#
# Env: TS_LIB, REPO (checkout dir, for the vendored fixture), OUTDIR,
#      NSEED_A (30), NSEED_B (15), BUDGETS_A ("90 180"), BUDGETS_B ("60 120"),
#      DATASETS_B, SHARD (SLURM_ARRAY_TASK_ID), NSHARD (SLURM_ARRAY_TASK_COUNT).
# ----------------------------------------------------------------------------
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("TS_LIB", .libPaths()[1]))
  library(TreeTools)
})

REPO    <- Sys.getenv("REPO", ".")
OUTDIR  <- Sys.getenv("OUTDIR", "dev/benchmarks")
NSEED_A <- as.integer(Sys.getenv("NSEED_A", "30"))
NSEED_B <- as.integer(Sys.getenv("NSEED_B", "15"))
BUDG_A  <- as.integer(strsplit(trimws(Sys.getenv("BUDGETS_A", "90 180")), "\\s+")[[1]])
BUDG_B  <- as.integer(strsplit(trimws(Sys.getenv("BUDGETS_B", "60 120")), "\\s+")[[1]])
DATA_B  <- strsplit(trimws(Sys.getenv("DATASETS_B",
             "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
SHARD   <- as.integer(Sys.getenv("SHARD", Sys.getenv("SLURM_ARRAY_TASK_ID", "0")))
NSHARD  <- as.integer(Sys.getenv("NSHARD", Sys.getenv("SLURM_ARRAY_TASK_COUNT", "1")))

# --- gaps-as-missing recode (self-contained copy of TreeSearch:::.GapsAsMissing) -
# The vendored fixture's two-island pathology only manifests under gaps-as-missing
# Fitch.  We recode at the DATA level (contrast-matrix: any gap-bearing token ->
# fully missing) and score with the default engine, which is byte-identical to the
# engine's own inapplicable = "missing" path (MaximizeParsimony.R: if "missing"
# then .GapsAsMissing(dataset); inapplicable <- "bgs").  Doing it here keeps the
# benchmark self-contained -- it does not depend on the built engine exposing the
# "missing" option (which post-dates the grafted branch this sweep builds from).
.GapsAsMissing <- function(dataset) {
  gapCol <- match("-", attr(dataset, "levels"))
  if (is.na(gapCol)) return(dataset)
  contrast <- attr(dataset, "contrast")
  contrast[contrast[, gapCol] == 1, ] <- 1
  attr(dataset, "contrast") <- contrast
  attr(dataset, "min.length") <- NULL
  dataset
}

# --- VARIANTS: dots overrides applied on top of strategy = "thorough" ---------
# base = current thorough; the rest probe the completeness/merge levers.
VARIANTS <- list(
  base        = list(),                                          # thorough as-is
  ws5         = list(wagnerStarts = 5L),                         # = intensive
  drift2      = list(driftCycles  = 2L),                         # uphill tunnelling
  drift5      = list(driftCycles  = 5L),
  ratchet10   = list(ratchetCycles = 10L),                       # cheaper reps -> throughput
  drift2_ws5  = list(driftCycles  = 2L, wagnerStarts = 5L)
)

# --- two-island fixture + golden references (vendored in the repo) -------------
fixt_dir <- file.path(REPO, "dev/benchmarks/data")
OG <- c("Galeaspida", "Osteostraci")
key1 <- function(t) ape::write.tree(SortTree(RootTree(t, OG)))
datA <- .GapsAsMissing(ReadAsPhyDat(file.path(fixt_dir, "zhu2013_2island.nex")))
k_isl2 <- vapply(ape::read.tree(file.path(fixt_dir, "zhu2013_2island_island2.tre")),
                 key1, character(1))
k_main <- vapply(ape::read.tree(file.path(fixt_dir, "zhu2013_2island_main.tre")),
                 key1, character(1))

# --- Part B datasets: classic benchmarks, gaps-as-missing Fitch ---------------
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
TARGET <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)

# --- enumerate all cells, then take this shard's stripe ------------------------
cellsA <- expand.grid(part = "A", variant = names(VARIANTS), budget = BUDG_A,
                      seed = seq_len(NSEED_A), dataset = "zhu2013_2island",
                      stringsAsFactors = FALSE)
cellsB <- expand.grid(part = "B", variant = names(VARIANTS), budget = BUDG_B,
                      seed = seq_len(NSEED_B), dataset = DATA_B,
                      stringsAsFactors = FALSE)
cells <- rbind(cellsA, cellsB)
mine  <- cells[(seq_len(nrow(cells)) - 1L) %% NSHARD == SHARD, , drop = FALSE]
cat(sprintf("shard %d/%d: %d of %d cells\n", SHARD, NSHARD, nrow(mine), nrow(cells)))

run_cell <- function(r) {
  dots <- VARIANTS[[r$variant]]
  set.seed(r$seed)
  if (r$part == "A") {
    # datA already recoded gaps-as-missing above; default engine = Fitch here.
    call <- c(list(datA, concavity = Inf,
                   strategy = "thorough", maxReplicates = 80L, targetHits = 9999L,
                   maxSeconds = r$budget, poolMaxSize = 2000L,
                   enumTimeFraction = 0.3, nThreads = 1L, verbosity = 0L), dots)
    t <- system.time(res <- suppressWarnings(do.call(MaximizeParsimony, call)))
    kr <- vapply(res, key1, character(1))
    data.frame(r, score = min(as.double(attr(res, "score"))),
               island2 = sum(k_isl2 %in% kr), main = sum(k_main %in% kr),
               both = sum(k_isl2 %in% kr) > 0 && sum(k_main %in% kr) > 0,
               n_topologies = length(kr),
               replicates = attr(res, "replicates"),
               over = NA_real_, elapsed = round(as.double(t["elapsed"]), 1))
  } else {
    phy <- fitch(inapplicable.phyData[[r$dataset]])
    call <- c(list(phy, strategy = "thorough", maxSeconds = r$budget,
                   nThreads = 1L, verbosity = 0L), dots)
    t <- system.time(res <- suppressWarnings(do.call(MaximizeParsimony, call)))
    sc <- min(as.double(attr(res, "score")))
    data.frame(r, score = sc, island2 = NA_integer_, main = NA_integer_,
               both = NA, n_topologies = length(res),
               replicates = attr(res, "replicates"),
               over = sc - TARGET[[r$dataset]],
               elapsed = round(as.double(t["elapsed"]), 1))
  }
}

# Accumulate results.  NB: pre-allocate and NEVER do `rows[[i]] <- NULL` -- in R
# that DELETES the element and shrinks the list, so the next `rows[[i]]` read
# runs off the end (subscript out of bounds) and kills the whole shard on the
# first errored cell.  We assign only successes; do.call(rbind, .) drops the
# untouched NULL slots.
rows <- vector("list", nrow(mine))
for (i in seq_len(nrow(mine))) {
  r <- mine[i, , drop = FALSE]
  rr <- tryCatch(run_cell(r), error = function(e) {
    cat(sprintf("  ERROR %s/%s b=%d s=%d: %s\n", r$part, r$variant, r$budget,
                r$seed, conditionMessage(e))); NULL })
  if (!is.null(rr)) {
    rows[[i]] <- rr
    cat(sprintf("  %s %-11s %-14s b=%3d s=%2d -> score=%.0f%s reps=%s [%.0fs]\n",
        rr$part, rr$variant, rr$dataset, rr$budget, rr$seed, rr$score,
        if (rr$part == "A") sprintf(" isl2=%d/16 both=%s", rr$island2, rr$both)
        else sprintf(" over=%+.0f", rr$over), rr$replicates, rr$elapsed))
  }
}
df <- do.call(rbind, rows)
out <- file.path(OUTDIR, sprintf("two_island_sweep_shard%03d.csv", SHARD))
if (is.null(df)) {
  cat(sprintf("shard %d: ALL %d cells failed; no CSV written\n", SHARD, nrow(mine)))
} else {
  write.csv(df, out, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", out, nrow(df)))
}
