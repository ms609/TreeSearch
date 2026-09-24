# D1 ORACLE (audit D1) — scoring-only, NO reinsertion, ZERO topology-risk.
# From an identical TNT `mult` T0 (global-TBR optimum), run our rss sectorial and
# let the C++ TS_FREE_HTU_PROBE diagnostic (ts_sector.cpp search_sector) compare,
# per sector: the HTU-ANCHORED frozen reduced score vs an UNCONSTRAINED reduced
# search where the HTU floats as an ordinary (S+1)th leaf (free re-resolve x
# re-attach). Since reduced = full - const (const invariant to HTU attachment),
# free < frozen on ANY sector PROVES a strictly shorter FULL tree exists that our
# anchored sectorial cannot reach -> D1 confirmed (the escape lever). If no sector
# shows free < frozen, D1 is refuted for the EW case.
#   advisor: the escape only shows via the from-scratch RAS rebuild w/ free HTU,
#   so rasStarts>=3 and the free probe does its own Wagner+TBR (R=3).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
Sys.setenv(TS_FREE_HTU_PROBE = "1")
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006")), "\\s+")[[1]]

get_t0 <- function(phy, seed = 1) {
  wd <- file.path(tempdir(), paste0("d1t0", Sys.getpid()))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "hold 100;", sprintf("rseed %d;", seed),
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
             file.path(wd, "dttest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  invisible(suppressWarnings(system2(TNT, args = "dttest.run;", stdout = TRUE, stderr = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  t0
}

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  t0 <- get_t0(phy)
  cat(sprintf("\n==== %s | T0 len=%.0f | D1 warm-revert probe (rss-only, rasStarts=1) ====\n",
              nm, TreeLength(t0, phy)))
  set.seed(1)
  invisible(suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
    nThreads = 1L, maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
    xssRounds = 0L, cssRounds = 0L, rssRounds = 1L, rasStarts = 20L, wagnerStarts = 1L,
    sectorMinSize = 30L, sectorMaxSize = 45L, fuseInterval = 9999L)))
}
