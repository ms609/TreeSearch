#!/usr/bin/env Rscript
# ============================================================================
# REGRESSION EXEMPLAR: "thorough" must capture BOTH MPT islands of Zhu2013.
# ============================================================================
#
# The original Zhu et al. (2013) Entelognathus matrix (Nature SI extraction,
# vendored here as data/zhu2013_2island.nex; NOT the cleaned-up bundled
# inapplicable.phyData[["Zhu2013"]], which no longer shows the pathology),
# scored with gaps-as-missing under Fitch parsimony, has its most-parsimonious
# trees (length 598) split across two TBR-plateau islands that are DISCONNECTED
# at equal score:
#
#   * a large "main" island that every search reaches readily, and
#   * a small island of exactly 16 trees ("island 2") whose basin is narrow:
#     seeding TreeSearch from any of the 16 and running a pure accept-equal TBR
#     walk closes over exactly those 16 and NONE of the main island (verified
#     2026-06-24). TS only reaches island 2 if a replicate's random start lands
#     in its basin. TNT (mult=replic 100 tbr; bbreak=tbr) finds it; a default
#     TreeSearch run frequently does NOT.
#
# A "thorough" search exists precisely to solve this completeness problem, so
# almost every run should recover BOTH islands. If a future "thorough" preset
# stops doing so, we have regressed inappropriately -- this benchmark fails.
#
# It is NOT a unit test (too slow, stochastic): run it on Hamilton, or locally
# with a small seed count as a smoke check.
#
# Usage:
#   Rscript dev/benchmarks/bench_thorough_two_islands.R            # full (Hamilton)
#   Rscript dev/benchmarks/bench_thorough_two_islands.R 2 60 20    # quick smoke
# Positional args / env overrides:
#   n_seeds        TS_BENCH_SEEDS       (default 10)
#   maxSeconds     TS_BENCH_SECONDS     (default 150)
#   maxReplicates  TS_BENCH_REPLICATES  (default 50)
# Other env:
#   TS_BENCH_THRESHOLD  fraction of seeds that must capture both islands (0.9)
#   TS_BENCH_STRATEGY   strategy preset under test (default "thorough")
# ----------------------------------------------------------------------------

suppressMessages({
  ok <- requireNamespace("TreeSearch", quietly = TRUE) &&
        !is.null(tryCatch(getNamespace("TreeSearch"), error = function(e) NULL))
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)        # dev checkout: use the source tree
  } else {
    library(TreeSearch)
  }
  library(TreeTools)
})

args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, env, default) {
  v <- if (length(args) >= i) args[[i]] else Sys.getenv(env, "")
  if (nzchar(v)) as.numeric(v) else default
}
N_SEEDS   <- as.integer(arg_or(1, "TS_BENCH_SEEDS", 10))
MAX_SECS  <- arg_or(2, "TS_BENCH_SECONDS", 150)
MAX_REPS  <- as.integer(arg_or(3, "TS_BENCH_REPLICATES", 50))
THRESHOLD <- as.numeric(Sys.getenv("TS_BENCH_THRESHOLD", "0.9"))
STRATEGY  <- Sys.getenv("TS_BENCH_STRATEGY", "thorough")

here   <- "dev/benchmarks/data"
og     <- c("Galeaspida", "Osteostraci")   # rooting outgroup for canonical keys
SCORE  <- 598                              # best parsimony length (gaps-as-missing)
dat    <- ReadAsPhyDat(file.path(here, "zhu2013_2island.nex"))
island2 <- ape::read.tree(file.path(here, "zhu2013_2island_island2.tre"))  # rare 16
main    <- ape::read.tree(file.path(here, "zhu2013_2island_main.tre"))     # main 256

# Canonical unrooted-topology key (rooting on the outgroup + SortTree makes it
# invariant to the search's arbitrary rooting / node order).
key1 <- function(t) ape::write.tree(SortTree(RootTree(t, og)))
k_island2 <- vapply(island2, key1, character(1))   # the rare 16
k_main    <- vapply(main,    key1, character(1))    # main island (256)

cat(sprintf(paste0("Two-island regression: strategy=%s | seeds=%d | %gs | ",
                   "maxReplicates=%d | island2=%d trees | threshold=%.0f%%\n"),
            STRATEGY, N_SEEDS, MAX_SECS, MAX_REPS, length(k_island2),
            100 * THRESHOLD))

rows <- vector("list", N_SEEDS)
for (s in seq_len(N_SEEDS)) {
  set.seed(s)
  res <- MaximizeParsimony(
    dat, concavity = Inf, inapplicable = "missing",
    strategy = STRATEGY, maxReplicates = MAX_REPS, targetHits = 9999L,
    maxSeconds = MAX_SECS,
    # poolMaxSize/enumTimeFraction overrides only -- the thorough preset still
    # governs the actual search heuristics. Pool must be able to HOLD both
    # islands, and enumeration must have time to expand a seeded island, so a
    # miss reflects SEEDING, not a capacity/time artefact.
    control = SearchControl(poolMaxSize = 2000L, enumTimeFraction = 0.3),
    verbosity = 0L, nThreads = 1L)

  kr <- vapply(res, key1, character(1))
  n2 <- sum(k_island2 %in% kr)           # how many of the 16 recovered
  nm <- sum(k_main    %in% kr)
  both <- n2 > 0 && nm > 0
  score <- attr(res, "score")
  rows[[s]] <- data.frame(seed = s, score = score, pool = length(kr),
                          island2 = n2, main = nm, both = both)
  cat(sprintf("  seed %2d: score=%.0f pool=%4d | island2 %2d/16 | main %3d | both=%s\n",
              s, score, length(kr), n2, nm, both))
}
df <- do.call(rbind, rows)
write.csv(df, "dev/benchmarks/results_two_islands.csv", row.names = FALSE)

frac <- mean(df$both)
full <- mean(df$island2 == length(k_island2))
pass <- frac >= THRESHOLD
cat(sprintf("\n==== %s ====\n", STRATEGY))
cat(sprintf("captured BOTH islands : %d/%d seeds (%.0f%%)  [target >= %.0f%%]\n",
            sum(df$both), N_SEEDS, 100 * frac, 100 * THRESHOLD))
cat(sprintf("fully recovered island2: %d/%d seeds (all 16)\n",
            sum(df$island2 == length(k_island2)), N_SEEDS))
cat(sprintf("RESULT: %s\n", if (pass) "PASS" else "FAIL (regression: thorough no longer reliably finds both islands)"))
quit(status = if (pass) 0L else 1L)
