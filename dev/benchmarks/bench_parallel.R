# PARALLEL BATCH runner — run a (dataset x seed) panel across a local PSOCK pool.
#
# For BATCH panels ONLY (e.g. an iterate-style panel over many seeds, or a
# preset sweep). NOT for a single authoritative candidate/timing measurement —
# oversubscription perturbs wall-clock and, under any wall-clock-bounded stop,
# the candidate count too. Each worker is single-threaded (nThreads=1, OMP=1)
# and REPLICATE-bounded, so candidates_evaluated stays valid per run.
#
# 8 physical cores, memory-bandwidth-bound Fitch -> realistic ~5-7x, not 16x.
# Set TS_HEADROOM high (>=4) while another panel/process is live.
#
# Env: TS_LIB, TS_DATASETS, TS_SEEDS, TS_REPS, TS_STRATEGY (auto), TS_HEADROOM, OUT_CSV.

suppressMessages(library(parallel))
LIB   <- normalizePath(Sys.getenv("TS_LIB", ".agent-p0"), winslash = "/")
WD    <- normalizePath(".", winslash = "/")
reps  <- as.integer(Sys.getenv("TS_REPS", "20"))
strat <- Sys.getenv("TS_STRATEGY", "auto")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5")), "\\s+")[[1]])
dsN   <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Wortley2006 Eklund2004 Zanol2014 Zhu2013 Giles2015 Dikow2009")), "\\s+")[[1]]
out   <- Sys.getenv("OUT_CSV", "dev/benchmarks/parallel_latest.csv")
headroom <- as.integer(Sys.getenv("TS_HEADROOM", "2"))
conc  <- max(1L, parallel::detectCores(logical = TRUE) - headroom)

jobs <- expand.grid(dataset = dsN, seed = seeds, stringsAsFactors = FALSE)
conc <- min(conc, nrow(jobs))
cat(sprintf("PARALLEL | %d jobs | conc=%d (cores=%d, headroom=%d) | %d reps | strategy=%s\n",
            nrow(jobs), conc, parallel::detectCores(logical = TRUE), headroom, reps, strat))

t0 <- Sys.time()
cl <- makePSOCKcluster(conc)
on.exit(stopCluster(cl))
clusterExport(cl, c("LIB", "WD", "reps", "strat", "jobs"), envir = environment())
invisible(clusterEvalQ(cl, {
  setwd(WD)                                    # PSOCK workers do NOT inherit CWD
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  suppressMessages({ library(TreeSearch, lib.loc = LIB); library(TreeTools) })
  data("inapplicable.phyData", package = "TreeSearch")
  fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
}))
rows <- parLapplyLB(cl, seq_len(nrow(jobs)), function(i) {
  nm <- jobs$dataset[i]; sd <- jobs$seed[i]
  d <- fitch(inapplicable.phyData[[nm]]); set.seed(sd)
  r <- suppressWarnings(MaximizeParsimony(d, maxReplicates = reps, targetHits = 999L,
                                          maxSeconds = 0, nThreads = 1L, strategy = strat,
                                          verbosity = 0L))
  data.frame(dataset = nm, seed = sd, score = attr(r, "score"),
             candidates = attr(r, "candidates_evaluated"), stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)
wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
agg <- do.call(rbind, lapply(split(res, res$dataset), function(d)
  data.frame(dataset = d$dataset[1], score_best = min(d$score),
             score_med = median(d$score), cand_med = median(d$candidates),
             stringsAsFactors = FALSE)))
cat(sprintf("done in %.0fs (%d jobs at conc=%d)\n", wall, nrow(jobs), conc))
print(agg[order(agg$dataset), ], row.names = FALSE)
write.csv(res, out, row.names = FALSE)
cat("rows ->", out, "\n")
