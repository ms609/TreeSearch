# FRAMING the TS-vs-TNT gap: efficiency x throughput decomposition.
#
# Answers the user's question BEFORE any VTune: "is TNT really running 10x
# more iter/sec, and is each iteration of similar complexity?"  Uses the
# bitness-independent candidate/rearrangement counts (local 32-bit TNT is
# fine for COUNTS) plus local wall to decompose the LOCAL wall ratio:
#
#   wall_ratio (TS/TNT) = efficiency (ts_cand / tnt_rearr)
#                       x throughput (tnt_rate / ts_rate)         [identity]
#
# CAVEAT (advisor): local TNT is 32-bit, so `throughput` is a LOWER BOUND on
# the true (64-bit) per-second advantage; ~2x correction noted separately.
# Counts are bitness-portable, so `efficiency` is exact.
#
# Run on a FRESH post-fix + unrooted-default build (.agent-p0).  All prior
# profiling data predates both the Wagner directional fix (2b299e4b) and the
# unrooted-TBR-default flip (25e35be7) and is stale.
#
# Env: TS_LIB (default .agent-p0), TNT_EXE, TS_SEEDS, TS_DATASETS, TS_SECONDS,
#      OUT_CSV.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})

TNT_EXE <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
secs    <- as.double(Sys.getenv("TS_SECONDS", "90"))
seeds   <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2")), "\\s+")[[1]])
replic  <- 50L
hits    <- 10L
out_csv <- Sys.getenv("OUT_CSV", "dev/profiling/framing_latest.csv")
dsNames <- strsplit(trimws(Sys.getenv("TS_DATASETS",
            "Wortley2006 Giles2015 Zhu2013 Zanol2014")), "\\s+")[[1]]

data("inapplicable.phyData", package = "TreeSearch")

tnt_work <- file.path(tempdir(), "tntframe")
dir.create(tnt_work, showWarnings = FALSE, recursive = TRUE)

run_tnt <- function(phy, seed, timeout_s) {
  datafile <- file.path(tnt_work, "datafile.tnt")
  runfile  <- file.path(tnt_work, "ftnt.run")
  logfile  <- file.path(tnt_work, "ftnt.log")
  WriteTntCharacters(phy, datafile)
  to <- sprintf("%02d:%02d:%02d", timeout_s %/% 3600,
                (timeout_s %% 3600) %/% 60, timeout_s %% 60)
  # `log <file>;` ... `log/;` captures TNT's screen output to a file.  Essential
  # on headless Linux: 64-bit TNT inits a curses UI (needs TERM set, e.g. xterm)
  # and uses the alternate-screen buffer, so the "Best score:" / "Total
  # rearrangements examined:" lines never reach a captured stdout pipe — but the
  # log file gets them regardless.  Works identically on Windows.
  script <- paste(
    "mxram 1024;",
    sprintf("log %s;", basename(logfile)),
    sprintf("proc %s;", basename(datafile)),
    "hold 10000;",
    sprintf("rseed %d;", seed),
    sprintf("timeout %s;", to),
    sprintf("xmult=hits %d replic %d;", hits, replic),
    "best;", "log/;", "quit;", sep = "\n")
  writeLines(script, runfile)
  old <- setwd(tnt_work); on.exit(setwd(old))
  if (file.exists(logfile)) file.remove(logfile)
  t0 <- Sys.time()
  # HEADLESS FIX (verified Hamilton 2026-06-19): feed the script via STDIN, not
  # as a run-file argument.  `tnt run.run;` makes 64-bit TNT try to EXECUTE the
  # filename as a command ("Must read data before generating random trees") and
  # `proc` never runs; piping commands to stdin is the canonical non-interactive
  # mode.  No curses is initialised when stdin is not a TTY, so TERM=dumb (set in
  # the sbatch) suffices — the earlier TERM=xterm / `log`-only theory was wrong.
  out <- tryCatch(
    system2(TNT_EXE, stdin = basename(runfile),
            stdout = TRUE, stderr = TRUE),
    error = function(e) character(0))
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  # Prefer the log file (robust under curses); fall back to captured stdout.
  txt <- if (file.exists(logfile)) readLines(logfile, warn = FALSE) else out
  txt <- iconv(txt, from = "", to = "UTF-8", sub = "")
  # This build prints "Best score (TBR): N" (not "Best score: N") — tolerate the
  # optional "(...)" between "Best score" and the colon.
  score <- suppressWarnings(as.double(
    sub(".*Best score[^:]*:\\s*([0-9.]+).*", "\\1",
        grep("Best score", txt, value = TRUE)[1])))
  rearr <- suppressWarnings(as.double(gsub(",", "",
    sub(".*Total rearrangements examined:\\s*([0-9,]+).*", "\\1",
        grep("Total rearrangements examined:", txt, value = TRUE)[1]))))
  list(score = score, rearr = rearr, wall = wall)
}

fitch_convert <- function(phy) {
  m <- PhyDatToMatrix(phy, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}

run_ts <- function(phy, seed, timeout_s) {
  set.seed(seed)
  t0 <- Sys.time()
  r <- suppressWarnings(MaximizeParsimony(
    phy, maxReplicates = 50L, nThreads = 1L, strategy = "auto",
    maxSeconds = timeout_s, verbosity = 0L))
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  list(score = attr(r, "score"),
       cand  = attr(r, "candidates_evaluated"),
       reps  = attr(r, "replicates"),
       wall  = wall)
}

cat(sprintf("FRAMING | %d datasets | seeds {%s} | cap %gs\n  TNT: %s\n",
            length(dsNames), paste(seeds, collapse = ","), secs, TNT_EXE))

rows <- list()
for (nm in dsNames) {
  raw   <- inapplicable.phyData[[nm]]
  fitch <- fitch_convert(raw)
  for (sd in seeds) {
    ts  <- run_ts(fitch, sd, secs)
    tnt <- run_tnt(fitch, sd, secs)
    eff  <- ts$cand / tnt$rearr
    tsr  <- ts$cand / ts$wall
    tntr <- tnt$rearr / tnt$wall
    thr  <- tntr / tsr
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, tips = length(raw), seed = sd,
      ts_score = ts$score, tnt_score = tnt$score,
      gapB = ts$score - tnt$score,
      ts_cand = ts$cand, tnt_rearr = tnt$rearr,
      efficiency = round(eff, 2),
      ts_rate = round(tsr / 1e6, 2), tnt_rate = round(tntr / 1e6, 2),
      throughput = round(thr, 2),
      ts_wall = round(ts$wall, 1), tnt_wall = round(tnt$wall, 1),
      wall_ratio = round(ts$wall / tnt$wall, 2),
      ts_reps = ts$reps, stringsAsFactors = FALSE)
    cat(sprintf("  %-12s s%d  TS %g(+%g) cand %.2gM  TNT %g rearr %.2gM | eff %.2f x thr %.2f = wall %.2f\n",
                nm, sd, ts$score, ts$score - tnt$score, ts$cand / 1e6,
                tnt$score, tnt$rearr / 1e6, eff, thr, ts$wall / tnt$wall))
  }
}
res <- do.call(rbind, rows)

agg <- do.call(rbind, lapply(split(res, res$dataset), function(d) data.frame(
  dataset = d$dataset[1], tips = d$tips[1],
  ts_best = min(d$ts_score), tnt_best = min(d$tnt_score),
  gapB = median(d$gapB),
  eff = median(d$efficiency),
  ts_rate_M = median(d$ts_rate), tnt_rate_M = median(d$tnt_rate),
  throughput = median(d$throughput),
  wall_ratio = median(d$wall_ratio),
  stringsAsFactors = FALSE)))
agg <- agg[order(agg$tips), ]

cat("\n=== Per-dataset decomposition (median; local 32-bit TNT) ===\n")
cat("eff = ts_cand/tnt_rearr (efficiency, bitness-free); throughput = tnt_rate/ts_rate (32-bit lower bound);\n")
cat("rates in M cand-or-rearr/sec; wall_ratio ~= eff x throughput (identity check).\n\n")
print(agg, row.names = FALSE)

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("\nPer-run rows -> %s\n", out_csv))
