# Discriminating trace for the TNT-audit acceptance-gate question:
# Does a sector ever improve on the REDUCED score (red_best < red_cur) while
# the FULL tree does NOT improve (full_new >= full_best)?
#
#   * If NEVER under EW-Fitch (- -> ?): reduced-improving <=> full-improving,
#     so the from-above HTU scoring is EXACT and the strict full-tree accept
#     gate is a NULL divergence from Goloboff 1999 (accept-on-reduced).
#   * If it HAPPENS under native NA (Brazeau, keep -): scoring is inexact and
#     the gate bites — but that is not the audited EW case.
#
# Uses the TS_SECT_DEBUG=1 trace already compiled into rss_search
# (ts_sector.cpp:1081). We force the sectorial path by giving a tree large
# enough for sectors and running with rss only.
#
# Env: TS_LIB (.agent-audit), TS_DS (dataset), TS_MODE (fitch|native), TS_SEED

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-audit"),
                                              winslash = "/"))
  library(TreeTools)
})

ds_name <- Sys.getenv("TS_DS", "Zanol2014")
mode    <- Sys.getenv("TS_MODE", "fitch")
seed    <- as.integer(Sys.getenv("TS_SEED", "1"))

data("inapplicable.phyData", package = "TreeSearch")
phy <- inapplicable.phyData[[ds_name]]
if (mode == "fitch") {
  m <- PhyDatToMatrix(phy, ambigNA = FALSE)
  m[m == "-"] <- "?"
  phy <- MatrixToPhyDat(m)
}

# Force sectorial to engage and run RSS specifically; disable ratchet/drift so
# the only score-changing phase whose trace we read is the sector accept path.
# rasStarts default 1; we test BOTH the polish (1) and rebuild (3) here.
ras <- as.integer(Sys.getenv("TS_RAS", "3"))
ctl <- SearchControl(
  ratchetCycles = 0L, driftCycles = 0L, nniPerturbCycles = 0L,
  pruneReinsertCycles = 0L, annealCycles = 0L,
  xssRounds = 0L, cssRounds = 0L,
  rssRounds = 3L, rasStarts = ras,
  sectorMinSize = 6L, sectorMaxSize = 50L,
  wagnerStarts = 1L, fuseInterval = 0L, intraFuse = FALSE,
  maxOuterResets = 0L, outerCycles = 1L
)

Sys.setenv(TS_SECT_DEBUG = "1")
set.seed(seed)
cat(sprintf("=== %s | mode=%s | seed=%d | rasStarts=%d ===\n",
            ds_name, mode, seed, ras))
# verbosity 0 so only the C-level TS_SECT_DEBUG REprintf lines appear on stderr
invisible(suppressWarnings(MaximizeParsimony(
  phy, maxReplicates = 1L, maxSeconds = 30, nThreads = 1L,
  strategy = "default", control = ctl, verbosity = 0L)))
cat("=== done ===\n")
