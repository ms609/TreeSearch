# Faithful plateau test (advisor): does letting the INTERNAL sector TBR hold many
# equal-length trees (sectorMaxHits) + keep laterals (sectorAcceptEqual) close the
# shared-start sectorial gap to TNT? TNT holds many trees while swapping a sector;
# we hold one (internal_max_hits = 1). rasStarts = 1 to isolate (re-solve is null).
#
# Engage-check: mh=20 must inflate candidates vs mh=1, else the knob isn't reaching
# the internal TBR. gap = TS_sect - TNT_sect from identical T0; lower = closer.
# STOP RULE: if 'plat' does not beat 'base', stop pulling levers -> instrument transfer.
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
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/sect_plateau.csv")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "sectplat"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

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
run_ts <- function(d, tree, rss, ras, aeq, mh) {
  set.seed(1)
  nt <- length(d)
  smin <- as.integer(round(nt * 0.35)); smax <- as.integer(round(nt * 0.65))
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, sectorMinSize = smin, sectorMaxSize = smax,
    rssRounds = as.integer(rss), rasStarts = as.integer(ras),
    sectorAcceptEqual = aeq, sectorMaxHits = as.integer(mh)))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s s%d: no T0\n", nm, sd)); next }
    base <- run_ts(phy, tn$t0, K, 1L, FALSE, 1L)    # current behaviour
    mh20 <- run_ts(phy, tn$t0, K, 1L, FALSE, 20L)   # +hold (max_hits alone)
    plat <- run_ts(phy, tn$t0, K, 1L, TRUE, 20L)    # +hold +accept_equal (faithful plateau)
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd,
      tnt = tn$s_sect, base = base$score, mh20 = mh20$score, plat = plat$score,
      g_base = base$score - tn$s_sect, g_mh20 = mh20$score - tn$s_sect,
      g_plat = plat$score - tn$s_sect,
      Mc_base = round(base$cand / 1e6, 1), Mc_plat = round(plat$cand / 1e6, 1),
      stringsAsFactors = FALSE)
    cat(sprintf("%-11s s%d | TNT=%.0f | base=%.0f mh20=%.0f plat=%.0f | gaps %+.0f/%+.0f/%+.0f | Mcand %.1f->%.1f\n",
                nm, sd, tn$s_sect, base$score, mh20$score, plat$score,
                base$score - tn$s_sect, mh20$score - tn$s_sect, plat$score - tn$s_sect,
                base$cand / 1e6, plat$cand / 1e6))
  }
}
S <- do.call(rbind, rows)
cat("\n== medians (gap = TS_sect - TNT_sect from identical T0; lower = closer) ==\n")
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], TNT = median(d$tnt),
  g_base = median(d$g_base), g_mh20 = median(d$g_mh20), g_plat = median(d$g_plat),
  Mc_base = median(d$Mc_base), Mc_plat = median(d$Mc_plat))))
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
