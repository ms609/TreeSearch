# Probe: SHARED-START sectorial quality — TNT vs TreeSearch from an IDENTICAL tree.
#
# Removes the starting-tree confound (the 2-31 step gap in bench_sectorial_yield.R)
# to isolate sectorial QUALITY. TNT builds ONE RAS+TBR tree T0 (hold 1), saves it
# parenthetically (tsave *), and runs sectsch=rss from it; we read T0 via
# ReadTntTree and run OUR sectorial from the SAME T0 (ratchet/drift OFF). Scores &
# rearrangement counts are bitness-independent, so the local 32-bit TNT is valid.
#
# sect_gap = TS_sect - TNT_sect, both from the same T0:
#   ~0  => our sectorial is quality-competitive (gap was the starting tree)
#   >0  => our sectorial is genuinely weaker (justifies the multi-start rewrite)
#
# Env: TS_LIB, TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS, OUT_CSV.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB",
            "dev/profiling/.vtune-lib-20260617081344"), winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
K   <- as.integer(Sys.getenv("TS_KPASS", "8"))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
         "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/sectorial_shared.csv")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "sectshared"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

run_tnt <- function(phy, seed, kpass) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", kpass), "quit;")
  writeLines(script, file.path(wd, "sharedstart.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "sharedstart.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_tbr <- num(sub(".*Best score \\(TBR\\):\\s*([0-9.]+).*", "\\1",
                   grep("Best score \\(TBR\\):", out, value = TRUE)[1]))
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  rearr <- num(sub(".*examined:\\s*([0-9,]+).*", "\\1",
                   grep("Total rearrangements examined:", out, value = TRUE)))
  t0 <- tryCatch(ReadTntTree(file.path(wd, "t0.tre")), error = function(e) NULL)
  if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, s_tbr = s_tbr,
       s_sect = if (length(s_sect)) s_sect[length(s_sect)] else NA,
       sect_rearr = if (length(rearr) >= 2) rearr[length(rearr)] - rearr[1] else NA)
}
run_ts <- function(d, tree, rss) {
  set.seed(1)
  nt <- length(d)
  smin <- as.integer(Sys.getenv("TS_SECTMIN", as.character(round(nt * 0.35))))
  smax <- as.integer(Sys.getenv("TS_SECTMAX", as.character(round(nt * 0.65))))
  use_css <- Sys.getenv("TS_USE_CSS", "0") == "1"
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L,
    cssRounds = if (use_css) as.integer(rss) else 0L,
    wagnerStarts = 1L, fuseInterval = 9999L,
    sectorMinSize = smin, sectorMaxSize = smax,
    rssRounds = if (use_css) 0L else as.integer(rss)))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s seed %d: no T0\n", nm, sd)); next }
    start <- TreeLength(tn$t0, phy)
    a <- run_ts(phy, tn$t0, 0L)   # TS TBR-only from T0 (control)
    b <- run_ts(phy, tn$t0, K)    # TS TBR + sectorial from T0
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd, start = start,
      tnt_sect = tn$s_sect, ts_tbr = a$score, ts_sect = b$score,
      sect_gap = b$score - tn$s_sect,
      tnt_Mrearr = round(tn$sect_rearr / 1e6, 2),
      ts_Mcand = round((b$cand - a$cand) / 1e6, 2), stringsAsFactors = FALSE)
  }
}
S <- do.call(rbind, rows)
cat("\n== Shared-start sectorial quality (TNT vs TS from identical T0) ==\n")
cat(sprintf("K=%d sectorial passes | seeds {%s}\n\n", K, paste(seeds, collapse = ",")))
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], start = median(d$start), TNT_sect = median(d$tnt_sect),
  TS_tbr = median(d$ts_tbr), TS_sect = median(d$ts_sect),
  sect_gap = median(d$sect_gap), stringsAsFactors = FALSE)))
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
