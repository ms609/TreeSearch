# SMOKE tier — breakage tripwire, ~seconds, run on every edit (POOL DRAINED).
#
# One R process, a few tiny datasets, REPLICATE-bounded (maxSeconds=0) so the
# candidate count is deterministic for a fixed seed (NOT wall-clock-bounded —
# that would make candidates machine-load-sensitive; see the critic note in
# dev/plans/2026-06-16-closing-the-tnt-gap.md). Green = "not broken / no
# candidate blow-up". This is a TRIPWIRE, never a ship gate: tiny datasets do
# not exercise sectorial search, so a real gap-lever can regress while smoke is
# green. Ship decisions use the iterate tier (bench_iterate.R).
#
# Env: TS_LIB (.agent-p0), TS_DATASETS, TS_REPS (4). SMOKE_WRITE_BASELINE=1 to
# (re)write dev/benchmarks/smoke_baseline.csv. Exit 1 on regression.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})
dsN  <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Longrich2010 Vinther2008 DeAssis2011")), "\\s+")[[1]]
reps <- as.integer(Sys.getenv("TS_REPS", "4"))
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }

t0 <- Sys.time()
res <- do.call(rbind, lapply(dsN, function(nm) {
  d <- fitch(inapplicable.phyData[[nm]]); set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(d, maxReplicates = reps, targetHits = 999L,
                                          maxSeconds = 0, nThreads = 1L, verbosity = 0L))
  data.frame(dataset = nm, score = attr(r, "score"),
             candidates = attr(r, "candidates_evaluated"), stringsAsFactors = FALSE)
}))
cat(sprintf("SMOKE | %d datasets | %d reps | %.1fs\n", length(dsN), reps,
            as.double(difftime(Sys.time(), t0, units = "secs"))))
print(res, row.names = FALSE)

base_f <- "dev/benchmarks/smoke_baseline.csv"
if (file.exists(base_f) && !nzchar(Sys.getenv("SMOKE_WRITE_BASELINE"))) {
  base <- read.csv(base_f)
  m <- merge(res, base, by = "dataset", suffixes = c("", ".base"))
  m$score_delta <- m$score - m$score.base
  m$cand_pct <- round(100 * (m$candidates / m$candidates.base - 1), 2)
  bad <- m[m$score_delta != 0 | abs(m$cand_pct) > 5, ]
  if (nrow(bad)) {
    cat("\nSMOKE FAIL (score changed or candidates moved >5%):\n")
    print(bad[, c("dataset", "score", "score.base", "cand_pct")], row.names = FALSE)
    quit(status = 1L)
  }
  cat("SMOKE OK (score unchanged; candidates within +/-5% of baseline)\n")
} else {
  write.csv(res, base_f, row.names = FALSE)
  cat("Wrote smoke baseline:", base_f, "\n")
}
