# CHEAP-WIN TEST: does ANCHORED rss at ~n/2 sector size + RAS restarts reach the
# sectsch target from T0 -- no floating/reinsert, pure selection+rasStarts tuning?
# ~n/2 is where free<anchored fired (D1 at the right size) AND anchored itself
# improved (563->560).  If anchored alone closes most of the gap -> cheap win
# (no kernel surgery); residual -> floating reinsert.  Shared start (TNT mult T0).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Wortley2006 Zhu2013 Giles2015")), "\\s+")[[1]]
target <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)

get_t0 <- function(phy) {
  wd <- file.path(tempdir(), paste0("nh", Sys.getpid()))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "hold 100;", "rseed 1;",
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
             file.path(wd, "nhtest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  invisible(suppressWarnings(system2(TNT, args = "nhtest.run;", stdout = TRUE, stderr = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  t0
}
rss_from <- function(phy, t0, ras, lo, hi) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L, nThreads = 1L,
        maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L,
        cssRounds = 0L, rssRounds = 8L, rasStarts = as.integer(ras), wagnerStarts = 1L,
        sectorMinSize = as.integer(lo), sectorMaxSize = as.integer(hi), fuseInterval = 9999L))
  min(as.double(attr(r, "score")))
}
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]]); n <- NTip(phy)
  lo <- max(6L, as.integer(round(n * 0.35))); hi <- as.integer(round(n * 0.55))
  t0 <- get_t0(phy); t0len <- TreeLength(t0, phy); tgt <- target[[nm]]
  cat(sprintf("\n==== %s (%dt, ~n/2 sector %d-%d) | T0=%.0f target=%d ====\n",
              nm, n, lo, hi, t0len, tgt))
  for (ras in c(3L, 10L, 20L)) {
    sc <- rss_from(phy, t0, ras, lo, hi)
    cat(sprintf("  anchored rss rasStarts=%-2d -> %.0f  (%+.0f vs T0, %+.0f vs target)%s\n",
                ras, sc, sc - t0len, sc - tgt, if (sc <= tgt) "  <== REACHED" else ""))
  }
}
