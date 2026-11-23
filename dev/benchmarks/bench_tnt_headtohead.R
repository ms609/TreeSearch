# TreeSearch vs TNT head-to-head — Phase 0 baseline harness.
#
# Establishes the authoritative, apples-to-apples picture against TNT 1.6:
#   * gap A (scoring method): TreeSearch Brazeau three-pass on RAW inapplicable
#     data vs TNT column-Fitch — NOT a search gap, shown for context.
#   * gap B (search quality):  TreeSearch Fitch (-> "?") vs TNT, same objective.
#   * candidates-per-improvement: TreeSearch `candidates_evaluated` (the new
#     instrumentation) vs TNT "Total rearrangements examined". Both quantities
#     are bitness-independent, so the LOCAL 32-bit TNT gives a valid comparison
#     of search efficiency; only wall-clock ratio needs 64-bit TNT (Hamilton).
#
# Modes (TS_MODE):
#   "converge" (default) — each engine runs to its natural convergence (capped
#       by TS_SECONDS as a safety timeout). Compares best score + candidate
#       counts. This is the candidates-per-improvement baseline.
#   "budget"             — both engines run to a fixed wall-clock (TS_SECONDS).
#       For the wall-clock ratio; only meaningful with a fair (64-bit) TNT.
#
# Env vars (all optional):
#   TS_LIB       library path for the instrumented TreeSearch build (.agent-p0)
#   TNT_EXE      path to tnt.exe (default: local 32-bit 1.6)
#   TS_DATASETS  space-separated dataset names from inapplicable.phyData
#   TS_SEEDS     space-separated integer seeds
#   TS_SECONDS   safety timeout (converge) or budget (budget mode), seconds
#   TS_MODE      "converge" | "budget"
#   TNT_REPLIC   TNT xmult replicates (default 50)
#   TNT_HITS     TNT xmult target hits  (default 10)
#   OUT_CSV      output CSV path (default dev/benchmarks/headtohead_latest.csv)

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})

TNT_EXE <- Sys.getenv("TNT_EXE",
                      "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
mode    <- Sys.getenv("TS_MODE", "converge")
secs    <- as.double(Sys.getenv("TS_SECONDS", "60"))
seeds   <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
replic  <- as.integer(Sys.getenv("TNT_REPLIC", "50"))
hits    <- as.integer(Sys.getenv("TNT_HITS", "10"))
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/headtohead_latest.csv")
dsNames <- strsplit(trimws(Sys.getenv("TS_DATASETS",
            "Wortley2006 Eklund2004 Zanol2014 Zhu2013 Giles2015 Dikow2009")),
            "\\s+")[[1]]

data("inapplicable.phyData", package = "TreeSearch")

# --- TNT helpers ----------------------------------------------------------
# TNT parses the script-name ARG as a command line (splits on digits/_), so the
# script basename must be purely alphabetic; data files read by `proc` are fine.
tnt_work <- file.path(tempdir(), "tntwork")
dir.create(tnt_work, showWarnings = FALSE, recursive = TRUE)

run_tnt <- function(phy, seed, timeout_s) {
  datafile <- file.path(tnt_work, "datafile.tnt")
  runfile  <- file.path(tnt_work, "htnt.run")
  WriteTntCharacters(phy, datafile)
  to <- sprintf("%02d:%02d:%02d", timeout_s %/% 3600,
                (timeout_s %% 3600) %/% 60, timeout_s %% 60)
  script <- paste(
    "mxram 1024;",
    sprintf("proc %s;", basename(datafile)),
    "hold 10000;",
    sprintf("rseed %d;", seed),
    sprintf("timeout %s;", to),
    sprintf("xmult=hits %d replic %d;", hits, replic),
    "best;", "quit;", sep = "\n")
  writeLines(script, runfile)
  old <- setwd(tnt_work); on.exit(setwd(old))
  t0 <- Sys.time()
  out <- tryCatch(
    system2(TNT_EXE, args = paste0(basename(runfile), ";"),
            stdout = TRUE, stderr = TRUE),
    error = function(e) character(0))
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  txt <- paste(out, collapse = "\n")
  score <- suppressWarnings(as.double(
    sub(".*Best score:\\s*([0-9.]+).*", "\\1",
        grep("Best score:", out, value = TRUE)[1])))
  rearr <- suppressWarnings(as.double(gsub(",", "",
    sub(".*Total rearrangements examined:\\s*([0-9,]+).*", "\\1",
        grep("Total rearrangements examined:", out, value = TRUE)[1]))))
  list(score = score, rearr = rearr, wall = wall)
}

# --- TreeSearch helper ----------------------------------------------------
fitch_convert <- function(phy) {
  m <- PhyDatToMatrix(phy, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}

run_ts <- function(phy, seed, timeout_s) {
  set.seed(seed)
  maxRep <- if (mode == "budget") 9999L else max(replic, 50L)
  t0 <- Sys.time()
  r <- suppressWarnings(MaximizeParsimony(
    phy, maxReplicates = maxRep, nThreads = 1L, strategy = "auto",
    maxSeconds = if (mode == "budget") timeout_s else timeout_s,
    verbosity = 0L))
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  list(score = attr(r, "score"),
       cand  = attr(r, "candidates_evaluated"),
       reps  = attr(r, "replicates"),
       wall  = wall)
}

# --- Run panel ------------------------------------------------------------
cat(sprintf("Head-to-head | mode=%s | %d datasets | seeds {%s} | cap %gs\n  TNT: %s\n",
            mode, length(dsNames), paste(seeds, collapse = ","), secs, TNT_EXE))
cat(strrep("-", 92), "\n")

rows <- list()
for (nm in dsNames) {
  raw   <- inapplicable.phyData[[nm]]
  fitch <- fitch_convert(raw)
  for (sd in seeds) {
    ts_f <- run_ts(fitch, sd, secs)
    ts_r <- run_ts(raw,   sd, secs)          # gap A: Brazeau three-pass
    tnt  <- run_tnt(fitch, sd, secs)
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, tips = length(raw), seed = sd,
      ts_fitch = ts_f$score, ts_raw = ts_r$score, tnt = tnt$score,
      gapB = ts_f$score - tnt$score,
      ts_cand = ts_f$cand, tnt_rearr = tnt$rearr,
      cand_ratio = round(ts_f$cand / tnt$rearr, 2),
      ts_wall = round(ts_f$wall, 1), tnt_wall = round(tnt$wall, 1),
      ts_reps = ts_f$reps, stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)

# --- Per-dataset summary --------------------------------------------------
agg <- do.call(rbind, lapply(split(res, res$dataset), function(d) {
  data.frame(
    dataset = d$dataset[1], tips = d$tips[1],
    ts_fitch_best = min(d$ts_fitch), ts_fitch_med = median(d$ts_fitch),
    tnt_best = min(d$tnt), gapB_med = median(d$gapB),
    ts_raw_med = median(d$ts_raw),                       # gap A context
    ts_cand_med = median(d$ts_cand), tnt_rearr_med = median(d$tnt_rearr),
    cand_ratio_med = median(d$cand_ratio),
    ts_wall_med = median(d$ts_wall), tnt_wall_med = median(d$tnt_wall),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(-agg$gapB_med), ]

print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("\nPer-run rows written to %s\n", out_csv))
