# Tier-2 A/B harness: deterministic same-seed runs of the standard-Fitch
# (TNT-parity) search. Score MUST be identical across builds (Tier 2 is
# behaviour-neutral); only wall-clock should move. nThreads=1, fixed seed.
#   env: TREESEARCH_VTUNE_LIB (lib path), TS_RUNS (default 10), TS_REPS (8)
LIBDIR <- Sys.getenv("TREESEARCH_VTUNE_LIB")
suppressMessages(library(TreeSearch, lib.loc = LIBDIR))
suppressMessages(library(TreeTools))

raw <- inapplicable.phyData[["Zhu2013"]]
m <- PhyDatToMatrix(raw, ambigNA = FALSE)
m[m == "-"] <- "?"                              # TNT-parity: standard Fitch
dataset <- MatrixToPhyDat(m)
stopifnot(!("-" %in% attr(dataset, "levels")))

N    <- as.integer(Sys.getenv("TS_RUNS", "10"))
reps <- as.integer(Sys.getenv("TS_REPS", "8"))

# One warm-up (page-in, allocator warm) excluded from stats.
set.seed(1)
invisible(suppressWarnings(MaximizeParsimony(dataset, maxReplicates = reps,
          nThreads = 1L, strategy = "auto", verbosity = 0L)))

times <- numeric(N); scores <- numeric(N)
for (i in seq_len(N)) {
  set.seed(1)                                   # identical work every run
  t0 <- proc.time()
  r <- suppressWarnings(MaximizeParsimony(dataset, maxReplicates = reps,
       nThreads = 1L, strategy = "auto", verbosity = 0L))
  times[i]  <- (proc.time() - t0)["elapsed"]
  scores[i] <- attr(r, "score")
}
cat(sprintf("LIB    : %s\n", LIBDIR))
cat(sprintf("scores : %s   (must be a single value)\n",
            paste(sort(unique(scores)), collapse = ",")))
cat(sprintf("median : %.3f s   (min %.3f / max %.3f / mean %.3f / sd %.3f)\n",
            median(times), min(times), max(times), mean(times), sd(times)))
cat(sprintf("all    : %s\n", paste(sprintf("%.2f", times), collapse = " ")))
