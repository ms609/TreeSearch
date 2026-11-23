# Phase-yield diagnosis (Phase 1) — where does TreeSearch spend its search?
#
# Uses existing instrumentation (no new build needed):
#   * attr(result, "timings")             per-phase cumulative wall-clock (ms)
#   * attr(result, "candidates_evaluated") total TBR/SPR candidates (Phase 0a)
#   * attr(result, "last_improved_rep")    replicate that last improved the best
#
# Localises the candidates-per-improvement gap to a phase BEFORE building any
# Phase 2 lever: which phase eats the wall-clock, and does the search keep
# improving late (effort well spent) or plateau early (effort wasted)?
#
# Sectorial = xss + rss + css. apples-to-apples Fitch (-> "?"), nThreads = 1.
#
# Env: TS_LIB (.agent-p0), TS_DATASETS, TS_SEEDS, TS_SECONDS (budget), OUT_CSV

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})

secs    <- as.double(Sys.getenv("TS_SECONDS", "30"))
seeds   <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/phase_yield_latest.csv")
dsNames <- strsplit(trimws(Sys.getenv("TS_DATASETS",
            "Wortley2006 Eklund2004 Zanol2014 Zhu2013 Giles2015 Dikow2009")),
            "\\s+")[[1]]

data("inapplicable.phyData", package = "TreeSearch")
fitch_convert <- function(phy) {
  m <- PhyDatToMatrix(phy, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}

rows <- list()
for (nm in dsNames) {
  fitch <- fitch_convert(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    set.seed(sd)
    r <- suppressWarnings(MaximizeParsimony(
      fitch, maxReplicates = 9999L, maxSeconds = secs, nThreads = 1L,
      strategy = "auto", verbosity = 0L))
    tm <- attr(r, "timings")
    g <- function(k) if (is.null(tm[[k]]) || is.na(tm[[k]])) 0 else tm[[k]]
    sect <- g("xss_ms") + g("rss_ms") + g("css_ms")
    total_ms <- sum(unlist(tm), na.rm = TRUE)
    reps <- attr(r, "replicates")
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, tips = length(fitch), seed = sd,
      score = attr(r, "score"),
      cand = attr(r, "candidates_evaluated"),
      reps = reps,
      last_improved = attr(r, "last_improved_rep"),
      # fraction of replicates AFTER the last improvement (= wasted effort)
      late_frac = round(1 - attr(r, "last_improved_rep") / max(reps, 1), 2),
      pct_wagner = round(100 * g("wagner_ms") / total_ms),
      pct_initial_tbr = round(100 * g("tbr_ms") / total_ms),
      pct_sector = round(100 * sect / total_ms),
      pct_ratchet = round(100 * g("ratchet_ms") / total_ms),
      pct_final_tbr = round(100 * g("final_tbr_ms") / total_ms),
      pct_fuse = round(100 * g("fuse_ms") / total_ms),
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)

agg <- do.call(rbind, lapply(split(res, res$dataset), function(d) {
  data.frame(
    dataset = d$dataset[1], tips = d$tips[1],
    score_med = median(d$score), cand_med = median(d$cand),
    reps_med = median(d$reps), late_frac_med = median(d$late_frac),
    sector = median(d$pct_sector), ratchet = median(d$pct_ratchet),
    final_tbr = median(d$pct_final_tbr), init_tbr = median(d$pct_initial_tbr),
    fuse = median(d$pct_fuse), wagner = median(d$pct_wagner),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(-agg$cand_med), ]
cat(sprintf("Phase-yield | %d datasets | seeds {%s} | %gs | nThreads=1\n",
            length(dsNames), paste(seeds, collapse = ","), secs))
cat("(phase columns = %% of wall-clock; late_frac = fraction of reps after last improvement)\n")
cat(strrep("-", 96), "\n")
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("\nPer-run rows written to %s\n", out_csv))
