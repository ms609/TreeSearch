# Phase 2 lever sweep — does cutting/rebalancing ratchet help the gap panel?
#
# Phase 1 found ratchet owns 63-83% of wall-clock (sectorial only 7-23%), the
# opposite of TNT. This tests ratchet/sectorial rebalancing via `auto` preset +
# `...` overrides (no rebuild). FIXED reps, parallel pool (replicate-bounded ->
# deterministic candidates; safe in the pool).
#
# CAVEAT: fixed-reps varies BOTH candidates and score per config, so it shows
# trade-offs, not a clean iso-candidate comparison (that needs the planned
# max_candidates C++ stop). Read: a config that holds score with FEWER candidates
# => ratchet over-invested; a config that improves score => quality win.
#
# Env: TS_LIB, TS_DATASETS, TS_SEEDS, TS_REPS, TS_HEADROOM, OUT_CSV.

suppressMessages(library(parallel))
LIB <- normalizePath(Sys.getenv("TS_LIB", ".agent-p0"), winslash = "/")
WD  <- normalizePath(".", winslash = "/")
reps  <- as.integer(Sys.getenv("TS_REPS", "20"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2")), "\\s+")[[1]])
dsN   <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Wortley2006 Eklund2004 Zanol2014 Zhu2013 Giles2015 Dikow2009")), "\\s+")[[1]]
headroom <- as.integer(Sys.getenv("TS_HEADROOM", "2"))

# Config sets, selectable via TS_SWEEP. Each value is a list of `...` overrides
# applied on top of strategy="auto". Round 1 (ratchet/sectorial) and round 2
# (fusing/ordering/starts) both gave no win over baseline — see
# dev/plans/2026-06-16-closing-the-tnt-gap.md Phase 2.
all_configs <- list(
  ratchet = list(
    baseline    = list(),
    ratchet6    = list(ratchetCycles = 6L),
    ratchet3    = list(ratchetCycles = 3L),
    adaptiveOff = list(adaptiveLevel = FALSE),
    sectorHeavy = list(xssRounds = 6L, rssRounds = 2L),
    rebalance   = list(ratchetCycles = 6L, xssRounds = 6L, rssRounds = 2L)
  ),
  fuse = list(
    baseline   = list(),
    intraFuse  = list(intraFuse = TRUE),
    fuseFreq   = list(fuseInterval = 1L),
    fuseEqual  = list(intraFuse = TRUE, fuseAcceptEqual = TRUE),
    clipTips   = list(clipOrder = 2L),
    wagner5    = list(wagnerStarts = 5L)
  ),
  optin = list(
    baseline  = list(),
    intraFuse = list(intraFuse = TRUE),
    wagner5   = list(wagnerStarts = 5L),
    combo     = list(intraFuse = TRUE, wagnerStarts = 5L)
  ),
  # Phase 3 probe: rebalance budget from ratchet toward EXACT sectorial (CSS),
  # which avoids the approximate XSS/RSS miss-and-revert waste. Tests whether
  # the cheapest exact phase, given more budget, carries more of the search.
  rebalance = list(
    baseline    = list(),
    css4        = list(cssRounds = 4L),
    ratchetDown = list(ratchetCycles = 8L),
    rebalA      = list(ratchetCycles = 12L, cssRounds = 4L),
    rebalB      = list(ratchetCycles = 8L, cssRounds = 4L, cssPartitions = 6L)
  )
)
configs <- all_configs[[Sys.getenv("TS_SWEEP", "ratchet")]]
if (is.null(configs)) stop("unknown TS_SWEEP")

jobs <- expand.grid(cfg = names(configs), dataset = dsN, seed = seeds,
                    stringsAsFactors = FALSE)
conc <- min(max(1L, parallel::detectCores(logical = TRUE) - headroom), nrow(jobs))
cat(sprintf("P2 levers | %d jobs (%d cfg x %d ds x %d seeds) | conc=%d | %d reps\n",
            nrow(jobs), length(configs), length(dsN), length(seeds), conc, reps))

t0 <- Sys.time()
cl <- makePSOCKcluster(conc)
on.exit(stopCluster(cl))
clusterExport(cl, c("LIB", "WD", "reps", "jobs", "configs"), envir = environment())
invisible(clusterEvalQ(cl, {
  setwd(WD); Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  suppressMessages({ library(TreeSearch, lib.loc = LIB); library(TreeTools) })
  data("inapplicable.phyData", package = "TreeSearch")
  fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
}))
rows <- parLapplyLB(cl, seq_len(nrow(jobs)), function(i) {
  cfg <- jobs$cfg[i]; nm <- jobs$dataset[i]; sd <- jobs$seed[i]
  d <- fitch(inapplicable.phyData[[nm]]); set.seed(sd)
  base <- list(d, strategy = "auto", maxReplicates = reps, targetHits = 999L,
               maxSeconds = 0, nThreads = 1L, verbosity = 0L)
  r <- suppressWarnings(do.call(MaximizeParsimony, c(base, configs[[cfg]])))
  data.frame(cfg = cfg, dataset = nm, seed = sd, score = attr(r, "score"),
             candidates = attr(r, "candidates_evaluated"), stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)
wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
agg <- aggregate(cbind(score, candidates) ~ cfg + dataset, res, median)
cat(sprintf("done in %.0fs\n", wall))
for (nm in dsN) {
  d <- agg[agg$dataset == nm, ]
  d <- d[order(d$score, d$candidates), ]
  b <- d[d$cfg == "baseline", ]
  d$dScore <- d$score - b$score
  d$dCand_pct <- round(100 * (d$candidates / b$candidates - 1))
  cat(sprintf("\n== %s (baseline %g @ %sM) ==\n", nm, b$score,
              format(round(b$candidates / 1e6), big.mark = ",")))
  print(d[, c("cfg", "score", "dScore", "dCand_pct")], row.names = FALSE)
}
write.csv(res, Sys.getenv("OUT_CSV", "dev/benchmarks/p2_levers.csv"), row.names = FALSE)
