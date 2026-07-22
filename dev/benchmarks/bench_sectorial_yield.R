# Probe 1: sectorial-search YIELD — TNT vs TreeSearch (measure the prize).
#
# Question: from a TBR-local-optimum start, how many steps does sectorial search
# buy, and at what rearrangement cost, in each engine? This sizes the prize for
# rewriting our sectorial (TNT constructs ~n/2 sectors + 3 RAS+TBR multi-start;
# we select [6,80] clades + single TBR — see dev/benchmarks/tnt_sector_defaults.csv).
#
# TNT:        `mult=replic 1` (RAS+TBR) -> K x `sectsch=rss`. Parses the per-phase
#             cumulative "Total rearrangements examined" + per-phase best score.
# TreeSearch: ONE replicate, ratchet/drift OFF, Wagner(1)+TBR; compare TBR-only
#             (rssRounds=0) vs TBR+RSS (rssRounds=K) on score + candidates_evaluated.
#
# Both rearrangement counts are bitness-independent, so the local 32-bit TNT is
# valid here (only wall-clock would need Hamilton). NB the two counters tally
# slightly different events (TNT counts all rearrangements incl. within-RAS; ours
# counts TBR/SPR candidates) so ABSOLUTE counts are only indicative — the
# steps-closed comparison is the clean signal.
#
# Env: TS_LIB, TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS, OUT_CSV.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB",
            "dev/profiling/.vtune-lib-20260617071429"), winslash = "/"))
  library(TreeTools)
})
TNT_EXE <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds   <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
K       <- as.integer(Sys.getenv("TS_KPASS", "5"))
dsN     <- strsplit(trimws(Sys.getenv("TS_DATASETS",
             "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/sectorial_yield.csv")

data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))

tnt_work <- file.path(tempdir(), "sectyield")
dir.create(tnt_work, showWarnings = FALSE, recursive = TRUE)

run_tnt_traj <- function(phy, seed, kpass) {
  datafile <- file.path(tnt_work, "datafile.tnt")
  runfile  <- file.path(tnt_work, "styield.run")
  WriteTntCharacters(phy, datafile)
  script <- c("mxram 1024;", sprintf("proc %s;", basename(datafile)),
              "hold 10000;", sprintf("rseed %d;", seed),
              "mult=replic 1;", rep("sectsch=rss;", kpass), "quit;")
  writeLines(script, runfile)
  old <- setwd(tnt_work); on.exit(setwd(old))
  out <- tryCatch(system2(TNT_EXE, args = paste0(basename(runfile), ";"),
                          stdout = TRUE, stderr = TRUE), error = function(e) character(0))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  rearr <- num(sub(".*examined:\\s*([0-9,]+).*", "\\1",
                   grep("Total rearrangements examined:", out, value = TRUE)))
  s_tbr <- num(sub(".*Best score \\(TBR\\):\\s*([0-9.]+).*", "\\1",
                   grep("Best score \\(TBR\\):", out, value = TRUE)[1]))
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  scores <- c(s_tbr, s_sect)
  n <- min(length(scores), length(rearr))
  if (n < 1) return(data.frame(phase=integer(0), score=double(0), cum_rearr=double(0)))
  data.frame(phase = 0:(n - 1), score = scores[1:n], cum_rearr = rearr[1:n])
}

run_ts <- function(phy, seed, kpass, do_sect) {
  set.seed(seed)
  r <- suppressWarnings(MaximizeParsimony(
    phy, maxReplicates = 1L, nThreads = 1L, strategy = "auto", maxSeconds = 0,
    verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
    xssRounds = 0L, cssRounds = 0L, wagnerStarts = 1L, fuseInterval = 9999L,
    rssRounds = if (do_sect) kpass else 0L))
  list(score = as.double(attr(r, "score")),
       cand  = as.double(attr(r, "candidates_evaluated")))
}

traj_all <- list(); summ <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tj <- run_tnt_traj(phy, sd, K); tj$dataset <- nm; tj$seed <- sd
    traj_all[[length(traj_all) + 1]] <- tj
    a <- run_ts(phy, sd, K, FALSE); b <- run_ts(phy, sd, K, TRUE)
    tnt_tbr  <- if (nrow(tj)) tj$score[1] else NA
    tnt_sect <- if (nrow(tj)) tj$score[nrow(tj)] else NA
    tnt_r0   <- if (nrow(tj)) tj$cum_rearr[1] else NA
    tnt_rF   <- if (nrow(tj)) tj$cum_rearr[nrow(tj)] else NA
    summ[[length(summ) + 1]] <- data.frame(
      dataset = nm, seed = sd,
      tnt_tbr = tnt_tbr, tnt_sect = tnt_sect, tnt_steps = tnt_tbr - tnt_sect,
      tnt_sect_Mrearr = round((tnt_rF - tnt_r0) / 1e6, 2),
      ts_tbr = a$score, ts_sect = b$score, ts_steps = a$score - b$score,
      ts_sect_Mcand = round((b$cand - a$cand) / 1e6, 2),
      stringsAsFactors = FALSE)
  }
}
traj <- do.call(rbind, traj_all); S <- do.call(rbind, summ)

cat("\n===== Sectorial yield: TNT vs TreeSearch (from TBR-local-opt) =====\n")
cat(sprintf("K=%d sectorial passes | seeds {%s}\n\n", K, paste(seeds, collapse = ",")))
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1],
  TNT_tbr = median(d$tnt_tbr), TNT_sect = median(d$tnt_sect),
  TNT_steps = median(d$tnt_steps), TNT_Mrearr = median(d$tnt_sect_Mrearr),
  TS_tbr = median(d$ts_tbr), TS_sect = median(d$ts_sect),
  TS_steps = median(d$ts_steps), TS_Mcand = median(d$ts_sect_Mcand),
  stringsAsFactors = FALSE)))
print(agg, row.names = FALSE)

cat("\n--- TNT sectorial score trajectory (median score by phase) ---\n")
for (nm in dsN) {
  t <- traj[traj$dataset == nm, ]
  if (!nrow(t)) next
  ph <- sort(unique(t$phase))
  med <- sapply(ph, function(p) median(t$score[t$phase == p]))
  cat(sprintf("  %-12s %s\n", nm, paste(sprintf("%g", med), collapse = " -> ")))
}

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
write.csv(traj, sub("\\.csv$", "_traj.csv", out_csv), row.names = FALSE)
cat(sprintf("\nWrote %s (+ _traj.csv)\n", out_csv))
