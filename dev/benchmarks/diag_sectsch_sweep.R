# SECTSCH LEVER SWEEP (parallel exploration of the OTHER sectsch differences).
# Shared-start: from TNT's mult T0, run our sectorial-ONLY (ratchet/drift OFF) under
# each lever config; report score reached vs the sectsch target.  Score-based +
# bounded (rssRounds fixed) => robust to CPU contention, deterministic.  Isolates
# which param-exposed sectsch lever (size, RAS restarts, accept-equal, max-hits)
# moves us toward TNT's sectsch endpoint, before any kernel work.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Wortley2006 Zhu2013 Giles2015")), "\\s+")[[1]]
target <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)  # TNT sectsch endpoints

get_t0 <- function(phy, seed = 1) {
  wd <- file.path(tempdir(), paste0("swt0", Sys.getpid(), substr(deparse(substitute(phy)),1,3)))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "hold 100;", sprintf("rseed %d;", seed),
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
             file.path(wd, "swtest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  invisible(suppressWarnings(system2(TNT, args = "swtest.run;", stdout = TRUE, stderr = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  t0
}
run_sect <- function(phy, t0, cfg) {
  set.seed(1)
  base <- list(dataset = phy, tree = t0, maxReplicates = 1L, nThreads = 1L,
               maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
               xssRounds = 0L, cssRounds = 0L, rssRounds = 8L, wagnerStarts = 1L,
               fuseInterval = 9999L)
  args <- modifyList(base, cfg)
  r <- tryCatch(suppressWarnings(do.call(MaximizeParsimony, args)),
                error = function(e) { message("ERR ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NA_real_)
  min(as.double(attr(r, "score")))
}
cfgs <- list(
  baseline    = list(),
  ras3        = list(rasStarts = 3L),
  ras6        = list(rasStarts = 6L),
  bigSectors  = list(sectorMinSize = 30L, sectorMaxSize = 45L),
  acceptEq    = list(sectorAcceptEqual = TRUE, sectorMaxHits = 10L),
  tntFaithful = list(rasStarts = 3L, sectorMinSize = 30L, sectorMaxSize = 45L,
                     sectorAcceptEqual = TRUE, sectorMaxHits = 10L)
)
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  t0 <- get_t0(phy); t0len <- TreeLength(t0, phy); tgt <- target[[nm]]
  cat(sprintf("\n==== %s | TNT mult T0=%.0f | sectsch target=%d ====\n", nm, t0len, tgt))
  for (cn in names(cfgs)) {
    sc <- run_sect(phy, t0, cfgs[[cn]])
    cat(sprintf("  %-12s -> %.0f  (%+.0f vs T0, %+.0f vs target)%s\n",
                cn, sc, sc - t0len, sc - tgt, if (is.finite(sc) && sc <= tgt) "  <== REACHED" else ""))
  }
}
