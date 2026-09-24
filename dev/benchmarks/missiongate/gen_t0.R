# Stage 1 (LOCAL, needs 32-bit TNT): generate the identical TNT `mult` T0 for the
# shared-T0 sectorial mission-gate, and record TNT's own sectsch-from-T0 target.
#
# For each (dataset, seed): TNT builds ONE RAS+TBR tree (hold 1, mult=replic 1),
# saves it (t0/<ds>_<seed>.tre), then runs sectsch=rss K times to record TNT's
# reachable sectorial score from that same T0. Scores are bitness-independent, so
# local 32-bit TNT is authoritative for the TARGET even though the sweep runs on
# Hamilton. Writes t0/manifest.csv (dataset, seed, start, tnt_tbr, tnt_sect).
#
# Mirrors bench_sectorial_shared.R's run_tnt(); EW-Fitch (- -> ?), no inapplicables.
#
# Env: TS_LIB (for TreeLength), TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS, OUT_DIR.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".missiongate-lib"),
                                              winslash = "/"))
  library(TreeTools)
})
TNT   <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
K     <- as.integer(Sys.getenv("TS_KPASS", "8"))
dsN   <- strsplit(trimws(Sys.getenv("TS_DATASETS",
           "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
outdir <- Sys.getenv("OUT_DIR", "dev/benchmarks/missiongate/t0")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num   <- function(x) as.double(gsub(",", "", x))
wd    <- file.path(tempdir(), "gent0"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

run_tnt <- function(phy, seed, kpass) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", kpass), "quit;")
  writeLines(script, file.path(wd, "sharedstart.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "sharedstart.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_tbr  <- num(sub(".*Best score \\(TBR\\):\\s*([0-9.]+).*", "\\1",
                    grep("Best score \\(TBR\\):", out, value = TRUE)[1]))
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- tryCatch(ReadTntTree(file.path(wd, "t0.tre")), error = function(e) NULL)
  if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, s_tbr = s_tbr,
       s_sect = if (length(s_sect)) s_sect[length(s_sect)] else NA)
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s seed %d: no T0\n", nm, sd)); next }
    start <- TreeLength(tn$t0, phy)
    # Save T0 with tip LABELS (species names), so Hamilton reads it label-safe.
    tf <- file.path(outdir, sprintf("%s_s%d.tre", nm, sd))
    ape::write.tree(tn$t0, tf)
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, seed = sd, nTip = length(phy), start = start,
      tnt_tbr = tn$s_tbr, tnt_sect = tn$s_sect,
      t0_file = basename(tf), stringsAsFactors = FALSE)
    cat(sprintf("%-12s s%d  start=%.0f  tnt_tbr=%.0f  tnt_sect=%.0f\n",
                nm, sd, start, tn$s_tbr, tn$s_sect))
  }
}
M <- do.call(rbind, rows)
write.csv(M, file.path(outdir, "manifest.csv"), row.names = FALSE)
cat(sprintf("\nWrote %d T0 trees + manifest to %s\n", nrow(M), outdir))
