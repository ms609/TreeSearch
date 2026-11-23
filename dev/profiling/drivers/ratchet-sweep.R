# Ratchet cost-effectiveness sweep. Ablation showed ratchet is load-bearing for
# quality (>=37t) but ~70-75% of wall. Two questions:
#   (a) Does FEWER ratchet cycles (12->6) hold quality at lower cost? (recipe tune)
#   (b) At MATCHED WALL, does ratchet-off + more replicates substitute for it?
# Configs vary ratchetCycles (named ... arg = override on auto preset) and
# maxReplicates. Report (wall, best score) per config => Pareto frontier.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Wills2012 Zanol2014 Zhu2013")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])

# label, ratchetCycles, maxReplicates
configs <- list(
  list(lab = "off_r6",  cyc = 0L,  reps = 6L),
  list(lab = "r6_r6",   cyc = 6L,  reps = 6L),
  list(lab = "r12_r6",  cyc = 12L, reps = 6L),   # current auto-ish
  list(lab = "off_r24", cyc = 0L,  reps = 24L)   # matched-wall substitution
)

trace("ts_driven_search", where = asNamespace("TreeSearch"),
      exit = quote(assign("TS_RES", list(timings = returnValue()$timings,
        replicates = returnValue()$replicates), envir = .GlobalEnv)),
      print = FALSE)

rows <- list()
for (nm in dsN) {
  phy <- fc(inapplicable.phyData[[nm]]); nt <- NTip(phy)
  for (sd in seeds) {
    for (cf in configs) {
      set.seed(sd); t0 <- Sys.time()
      r <- suppressWarnings(MaximizeParsimony(
        phy, maxReplicates = cf$reps, nThreads = 1L, strategy = "auto",
        verbosity = 0L, ratchetCycles = cf$cyc))
      w <- as.double(difftime(Sys.time(), t0, units = "secs"))
      res <- get("TS_RES", envir = .GlobalEnv); tm <- unlist(res$timings)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = nm, ntip = nt, seed = sd, cfg = cf$lab,
        score = attr(r, "score"), wall = round(w, 2),
        ratchet_ms = round(as.numeric(tm["ratchet_ms"])),
        sect_ms = round(sum(tm[c("xss_ms", "rss_ms", "css_ms")], na.rm = TRUE)),
        reps = res$replicates)
      x <- rows[[length(rows)]]
      cat(sprintf("%-12s nt=%d seed%d %-8s score=%s wall=%6.2fs reps=%s\n",
                  nm, nt, sd, cf$lab, x$score, x$wall, x$reps)); flush.console()
    }
  }
}
d <- do.call(rbind, rows)
write.csv(d, "dev/profiling/ratchet_sweep.csv", row.names = FALSE)

cat("\n=== Per dataset: best-known score, then per-config (misses / mean wall) ===\n")
for (nm in dsN) {
  sub <- d[d$dataset == nm, ]; best <- min(sub$score)
  cat(sprintf("\n%-12s best=%g\n", nm, best))
  for (cf in configs) {
    cc <- sub[sub$cfg == cf$lab, ]
    cat(sprintf("  %-8s scores %-18s miss %d/%d  mean wall %.2fs\n",
                cf$lab, paste(cc$score, collapse = ","),
                sum(cc$score > best), nrow(cc), mean(cc$wall)))
  }
}
cat("\nSWEEP COMPLETE\n")
