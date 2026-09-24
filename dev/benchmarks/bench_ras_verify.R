# Verify (with the now-plumbed rasStarts knob, clean .agent-ratchet build):
# does sector RE-SOLVE close the SHARED-START sectorial gap to TNT?
#
# Same design as bench_sectorial_shared.R: TNT builds ONE RAS+TBR tree T0 (hold 1),
# runs sectsch=rss from it; we read T0 and run OUR sectorial from the SAME T0 with
# ratchet/drift OFF, at rasStarts = 1 (polish) and rasStarts = 3 (re-solve). Scores
# are bitness-independent so local 32-bit TNT is valid. gap = TS_sect - TNT_sect.
#   re-solve gap ~ polish gap  => re-solve is NOT the missing piece (fidelity is)
#   re-solve gap << polish gap => re-solve closes it
#
# Env: TS_LIB (default .agent-ratchet), TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS, OUT_CSV.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
K   <- as.integer(Sys.getenv("TS_KPASS", "8"))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
         "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/ras_verify.csv")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "rasverify"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

run_tnt <- function(phy, seed, kpass) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", kpass), "quit;")
  writeLines(script, file.path(wd, "ss.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "ss.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- tryCatch(ReadTntTree(file.path(wd, "t0.tre")), error = function(e) NULL)
  if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, s_sect = if (length(s_sect)) s_sect[length(s_sect)] else NA)
}
run_ts <- function(d, tree, rss, ras) {
  set.seed(1)
  nt <- length(d)
  smin <- as.integer(round(nt * 0.35)); smax <- as.integer(round(nt * 0.65))
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, sectorMinSize = smin, sectorMaxSize = smax,
    rssRounds = as.integer(rss), rasStarts = as.integer(ras)))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s s%d: no T0\n", nm, sd)); next }
    start <- TreeLength(tn$t0, phy)
    a <- run_ts(phy, tn$t0, K, 1L)   # polish    (rasStarts = 1)
    b <- run_ts(phy, tn$t0, K, 3L)   # re-solve  (rasStarts = 3)
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd, start = start,
      tnt_sect = tn$s_sect, ts_polish = a$score, ts_resolve = b$score,
      gap_polish = a$score - tn$s_sect, gap_resolve = b$score - tn$s_sect,
      Mcand_polish = round(a$cand / 1e6, 2), Mcand_resolve = round(b$cand / 1e6, 2),
      stringsAsFactors = FALSE)
    cat(sprintf("%-11s s%d | start=%.0f TNT=%.0f | polish=%.0f(g%+.0f) resolve=%.0f(g%+.0f) | Mcand %.1f->%.1f\n",
                nm, sd, start, tn$s_sect, a$score, a$score - tn$s_sect,
                b$score, b$score - tn$s_sect, a$cand / 1e6, b$cand / 1e6))
  }
}
S <- do.call(rbind, rows)
cat("\n== medians (gap = TS_sect - TNT_sect, from identical T0) ==\n")
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], TNT = median(d$tnt_sect),
  polish = median(d$ts_polish), resolve = median(d$ts_resolve),
  gap_polish = median(d$gap_polish), gap_resolve = median(d$gap_resolve))))
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
