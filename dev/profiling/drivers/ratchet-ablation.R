# Ratchet-redundancy ablation: post-fix, does sectorial+TBR reach the same
# optimum WITHOUT the ratchet phase (which is ~60% of wall)?  Compares the
# `auto` recipe with ratchet ON vs ratchet disabled (ratchetCycles=0L passed as
# a named ... arg, which overrides ONLY that field on top of the auto preset --
# control=SearchControl(...) would wipe the preset, see MaximizeParsimony.R:593).
#
# Reports per (dataset, seed): score + wall + candidates for both configs, plus
# ratchet_ms / sectorial_ms (via trace on ts_driven_search) to confirm the
# disable bites.  Quality test = does ratchet-OFF match the best-known score?
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Zhu2013 Zanol2014")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
reps <- as.integer(Sys.getenv("TS_REPS", "8"))

trace("ts_driven_search", where = asNamespace("TreeSearch"),
      exit = quote(assign("TS_RES", list(
        timings = returnValue()$timings,
        replicates = returnValue()$replicates), envir = .GlobalEnv)),
      print = FALSE)

run1 <- function(phy, sd, ratchet_off) {
  set.seed(sd); t0 <- Sys.time()
  extra <- if (ratchet_off) list(ratchetCycles = 0L) else list()
  r <- suppressWarnings(do.call(MaximizeParsimony, c(list(
    phy, maxReplicates = reps, nThreads = 1L, strategy = "auto",
    verbosity = 0L), extra)))
  w <- as.double(difftime(Sys.time(), t0, units = "secs"))
  res <- get("TS_RES", envir = .GlobalEnv)
  tm <- unlist(res$timings)
  data.frame(score = attr(r, "score"), cand = attr(r, "candidates_evaluated"),
             wall = round(w, 2),
             ratchet_ms = round(as.numeric(tm["ratchet_ms"])),
             sect_ms = round(sum(tm[c("xss_ms", "rss_ms", "css_ms")], na.rm = TRUE)),
             reps = res$replicates)
}

rows <- list()
for (nm in dsN) {
  phy <- fc(inapplicable.phyData[[nm]])
  nt <- NTip(phy)
  for (sd in seeds) {
    for (off in c(FALSE, TRUE)) {
      x <- run1(phy, sd, off)
      x$dataset <- nm; x$ntip <- nt; x$seed <- sd
      x$cfg <- if (off) "off" else "on"
      rows[[length(rows) + 1]] <- x
      cat(sprintf("%-12s nt=%d seed%d ratchet=%-3s score=%s wall=%6.2fs ratchet_ms=%5s sect_ms=%5s reps=%s\n",
                  nm, nt, sd, x$cfg, x$score, x$wall, x$ratchet_ms, x$sect_ms, x$reps)); flush.console()
    }
  }
}
d <- do.call(rbind, rows)
write.csv(d, "dev/profiling/ratchet_ablation.csv", row.names = FALSE)

cat("\n=== Quality: does ratchet-OFF match best-known score per dataset? ===\n")
for (nm in dsN) {
  sub <- d[d$dataset == nm, ]
  best <- min(sub$score)
  on  <- sub[sub$cfg == "on", ]
  off <- sub[sub$cfg == "off", ]
  cat(sprintf("%-12s best=%g | ON  scores %s (miss %d/%d) | OFF scores %s (miss %d/%d) | wall ON=%.2fs OFF=%.2fs ratio=%.2f\n",
              nm, best,
              paste(on$score, collapse = ","), sum(on$score > best), nrow(on),
              paste(off$score, collapse = ","), sum(off$score > best), nrow(off),
              sum(on$wall), sum(off$wall), sum(off$wall) / sum(on$wall)))
}
cat("\nABLATION COMPLETE\n")
