# SECTORIAL LEG, cheap precursor (D2): does our FROZEN rss from TNT's T0 reach the
# sectsch target if we just raise rasStarts (TNT does R=3 + r=3)?  rss-only, from
# the identical TNT mult T0, bounded work (rssRounds fixed), nThreads=1.
#   reaches target  -> D2 (cheap, no kernel surgery)
#   stuck at ~T0    -> frozen rebuild is null; D1 (floating HTU) is the real lever
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014")), "\\s+")[[1]]
target <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)

get_t0 <- function(phy, seed = 1) {
  wd <- file.path(tempdir(), paste0("rrt0", Sys.getpid()))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "hold 100;", sprintf("rseed %d;", seed),
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
             file.path(wd, "rttest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  invisible(suppressWarnings(system2(TNT, args = "rttest.run;", stdout = TRUE, stderr = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  t0
}
rss_from <- function(phy, t0, ras) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
        nThreads = 1L, maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L,
        driftCycles = 0L, xssRounds = 0L, cssRounds = 0L, rssRounds = 8L,
        rasStarts = as.integer(ras), wagnerStarts = 1L, fuseInterval = 9999L))
  min(as.double(attr(r, "score")))
}
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  t0 <- get_t0(phy); t0len <- TreeLength(t0, phy); tgt <- target[[nm]]
  cat(sprintf("\n==== %s | TNT mult T0=%.0f | sectsch target=%d ====\n", nm, t0len, tgt))
  for (ras in c(1L, 3L, 6L)) {
    sc <- rss_from(phy, t0, ras)
    cat(sprintf("  rss rasStarts=%d -> %.0f  (%+.0f vs T0, %+.0f vs target)%s\n",
                ras, sc, sc - t0len, sc - tgt, if (sc <= tgt) "  <== REACHED" else ""))
  }
}
