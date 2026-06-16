# Standard-Fitch TNT-parity profiling driver — Area #5
# GOAL: profile the *standard Fitch* path that the TNT-parity benchmark uses.
#   TNT-parity replaces inapplicable "-" with missing "?" so both engines
#   optimise the identical Fitch objective (no Brazeau-Guillerme-Smith NA
#   handling).  Removing the "-" level makes the C++ engine take has_na=FALSE
#   and use the flat / 4-wide (T-245) kernels — a code path NEVER profiled
#   before (all prior rounds used the NA path on raw inapplicable.phyData).
#
# Reports per-phase timings via attr(result, "timings") so we can see where
# the standard-Fitch wall-clock actually goes, vs the (stale, NA-path) round-6
# distribution in dev/expertise/profiling.md.
#
# Params (env): TS_DATASET (default Zhu2013), TS_REPS (default 3), TS_SEED (1)
# nThreads=1 always (apples-to-apples with single-threaded TNT xmult).

LIBDIR <- Sys.getenv("TREESEARCH_VTUNE_LIB",
                     unset = "dev/profiling/.vtune-lib-20260616051420")
suppressMessages(library(TreeSearch, lib.loc = LIBDIR))
suppressMessages(library(TreeTools))

ds_name <- Sys.getenv("TS_DATASET", unset = "Zhu2013")
n_reps  <- as.integer(Sys.getenv("TS_REPS", unset = "3"))
seed    <- as.integer(Sys.getenv("TS_SEED", unset = "1"))

raw <- inapplicable.phyData[[ds_name]]

# --- Convert inapplicable "-" -> missing "?" (TNT-parity Fitch objective) ---
m <- PhyDatToMatrix(raw, ambigNA = FALSE)
n_dash <- sum(m == "-")
m[m == "-"] <- "?"
dataset <- MatrixToPhyDat(m)
lv <- attr(dataset, "levels")
stopifnot(!("-" %in% lv))  # confirm standard-Fitch path (has_na = FALSE)

cat(sprintf("Dataset: %s | %d tips | %d patterns | %d levels (%s) | %d '-' -> '?'\n",
            ds_name, length(dataset), attr(dataset, "nr"),
            length(lv), paste(lv, collapse = ""), n_dash))

# Auto strategy will pick a preset from nTip/nChar; report it.
strat <- TreeSearch:::.AutoStrategy(length(dataset), attr(dataset, "nr"))
cat(sprintf("Auto strategy -> %s\n", strat))

set.seed(seed)
t0 <- proc.time()
result <- suppressWarnings(
  MaximizeParsimony(
    dataset,
    maxReplicates = n_reps,
    nThreads      = 1L,
    strategy      = "auto",
    verbosity     = 0L
  )
)
elapsed <- (proc.time() - t0)["elapsed"]

tm <- attr(result, "timings")
tm <- tm[order(-tm)]
tot <- sum(tm, na.rm = TRUE)

cat(sprintf("\nElapsed: %.2f s | Score: %s | Reps: %s | MPTs: %s\n",
            elapsed, attr(result, "score"), attr(result, "replicates"),
            attr(result, "n_topologies")))
cat(sprintf("Sum of phase timings: %.1f ms\n\n", tot))
cat("Phase distribution (cumulative ms across all replicates):\n")
for (nm in names(tm)) {
  cat(sprintf("  %-22s %8.1f ms  %5.1f%%\n", nm, tm[[nm]],
              100 * tm[[nm]] / tot))
}
