# Isolated sectorial (RSS) profiling driver — Area #3
#
# Isolates ONLY the sector search component via the thin ts_rss_search Rcpp
# wrapper (no ratchet/fuse wrapper). Per the component-isolation plan
# (dev/plans/2026-06-19-component-isolation-profiling.md) and the advisor:
#   * Crank rssPicks HIGH and use FEW calls so the trailing global TBR and the
#     per-call make_dataset/init_from_edge marshaling (a DRIVER ARTIFACT) are
#     amortised — otherwise the profile is mostly TBR, not sectorial.
#   * EW-Fitch only ('-' -> '?'); the NA path is owned by another workstream.
#
# Start tree: a Wagner addition tree (non-optimal) so many sectors improve and
# the reinsert + full-tree rescore accept-path is exercised (bucket-2 coverage).
# Set TS_START=tbr to TBR-converge the start first (fewer accepts).
#
# Params (env): TS_DATASET (Zanol2014), TS_PICKS (80), TS_CALLS (12),
#   TS_SEED (1), TS_MINSIZE (6), TS_MAXSIZE (50), TS_START (wagner|tbr),
#   TS_ACCEPTEQ (0), TS_RATCHET (6 — internal_ratchet_cycles, unused by rss).
#
# bare target: <= 5 s. nThreads=1 (serial; sectorial RNG pulls from R's stream).

LIBDIR <- Sys.getenv("TREESEARCH_VTUNE_LIB", unset = ".agent-sect")
suppressMessages(library(TreeSearch, lib.loc = LIBDIR))
suppressMessages(library(TreeTools))

ds_name  <- Sys.getenv("TS_DATASET", unset = "Zanol2014")
n_picks  <- as.integer(Sys.getenv("TS_PICKS",   unset = "80"))
n_calls  <- as.integer(Sys.getenv("TS_CALLS",   unset = "12"))
seed     <- as.integer(Sys.getenv("TS_SEED",    unset = "1"))
min_size <- as.integer(Sys.getenv("TS_MINSIZE", unset = "6"))
max_size <- as.integer(Sys.getenv("TS_MAXSIZE", unset = "50"))
start_kind <- Sys.getenv("TS_START", unset = "wagner")
accept_eq <- as.integer(Sys.getenv("TS_ACCEPTEQ", unset = "0")) != 0L

raw <- inapplicable.phyData[[ds_name]]

# --- EW standard-Fitch objective: inapplicable '-' -> missing '?' ---
m <- PhyDatToMatrix(raw, ambigNA = FALSE)
m[m == "-"] <- "?"
dataset <- MatrixToPhyDat(m)
at <- attributes(dataset)
ds <- list(
  contrast = at$contrast,
  tip_data = matrix(unlist(dataset, use.names = FALSE),
                    nrow = length(dataset), byrow = TRUE),
  weight = at$weight,
  levels = at$levels
)
n_tip <- length(dataset)

# --- Build the start tree once (deterministic) ---
# Seed BEFORE the start build so the (RNG-using) ts_tbr_search start is
# reproducible across processes — required for the byte-identical A/B gate.
set.seed(seed)
wag <- TreeSearch:::ts_wagner_tree(ds$contrast, ds$tip_data, ds$weight, ds$levels)
start_edge <- wag$edge
start_score <- wag$score
if (identical(start_kind, "tbr")) {
  tb <- TreeSearch:::ts_tbr_search(start_edge, ds$contrast, ds$tip_data,
                                   ds$weight, ds$levels, maxHits = 1L)
  start_edge <- tb$edge
  start_score <- tb$score
}

cat(sprintf("Dataset: %s | %d tips | %d patterns | start(%s)=%g | picks=%d calls=%d\n",
            ds_name, n_tip, attr(dataset, "nr"), start_kind, start_score,
            n_picks, n_calls))

set.seed(seed)
scores <- numeric(n_calls)
n_searched <- integer(n_calls)
n_improved <- integer(n_calls)
t0 <- proc.time()
for (i in seq_len(n_calls)) {
  res <- TreeSearch:::ts_rss_search(
    start_edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    minSectorSize = min_size, maxSectorSize = max_size,
    acceptEqual = accept_eq, rssPicks = n_picks,
    ratchetCycles = 0L, maxHits = 1L)
  scores[i]     <- res$score
  n_searched[i] <- res$n_sectors_searched
  n_improved[i] <- res$n_sectors_improved
}
elapsed <- (proc.time() - t0)["elapsed"]

cat(sprintf("Elapsed: %.2f s | %d calls x %d picks | sectors searched=%d improved=%d\n",
            elapsed, n_calls, n_picks, sum(n_searched), sum(n_improved)))
cat(sprintf("Scores: min=%g median=%g max=%g\n",
            min(scores), median(scores), max(scores)))
# Gate signal (bit-identical A/B): per-call score + sector counts.
cat("GATE", paste(scores, collapse = ","), "|",
    paste(n_searched, collapse = ","), "|",
    paste(n_improved, collapse = ","), "\n")
