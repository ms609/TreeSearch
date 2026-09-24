# (c)-MECHANISM TEST (Goloboff sectsch escape), ZERO C++ risk -- all knobs R-exposed.
# Subagent (dev/plans/2026-06-17-sectsch-escape-mechanism.md) overturned D1 (frozen
# HTU; d1_confirm.out shows 0 confirms) and pins the escape on TNT's selectem GEOMETRY:
# LARGE (~n/2) sectors whose deep sub-clades are COLLAPSED into composite terminals, so
# RAS+TBR reshuffles whole sub-clades across the backbone -- large-radius full-tree moves
# our small-clade ras1 sectorial never proposes.  All this is reachable from R:
#   sectorMinSize/Max (force large) + sectorCollapseTarget (>0 collapse) + rasStarts(=3
#   RAS rebuild) + sectorAcceptEqual (the +1 bridge).  Walk-up-from-random selection is
#   the ONE piece NOT R-reachable (we pick existing in-band clades, not walk-up clades).
#
# LADDER (isolates each factor), SHARED start (TNT mult T0), 2 seeds, all-else-off:
#   base       defaults (min6 max50 ras1 coll0 eq F)        -- current behaviour
#   bigNoColl  min31 max99 ras3 coll0  eq F                 -- large sector, NO collapse
#   coll30     min31 max99 ras3 coll30 eq F                 -- + collapse to ~30 skeleton
#   coll30eq   min31 max99 ras3 coll30 eq T                 -- + accept-equal bridge
# collapse FIRES iff a picked clade > collapse_target; min31 > coll30 => every eligible
# pick collapses, so DIVERSITY (eligible clades in [31,99]) > 0 PROVES collapse fires
# (advisor's firing check, no rebuild).  Low diversity (1-2) + null => walk-up is the
# missing piece (implement next), NOT "(c) refuted" (pre-committed interpretation).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband2"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) suppressWarnings(as.double(gsub(",", "", x)))
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Wortley2006 Zhu2013 Giles2015")), "\\s+")[[1]]
target <- c(Wortley2006 = 479, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
ROUNDS <- as.integer(Sys.getenv("TS_RSSROUNDS", "15"))
SEEDS  <- as.integer(strsplit(Sys.getenv("TS_SEEDS", "1 2"), "\\s+")[[1]])

# rooted clade sizes (proxy for the C++ eligible set: subtree_size per internal node)
cladeSizes <- function(tree) {
  nTip <- length(tree$tip.label)
  po <- Postorder(tree)$edge
  cnt <- integer(max(po))
  cnt[seq_len(nTip)] <- 1L
  for (i in seq_len(nrow(po))) cnt[po[i, 1]] <- cnt[po[i, 1]] + cnt[po[i, 2]]
  cnt[(nTip + 1):length(cnt)]
}

# config = list(min, max, ras, coll, eq)
cfgs <- list(
  base      = list(6L,  50L, 1L,  0L, FALSE),
  bigNoColl = list(31L, 99L, 3L,  0L, FALSE),
  coll30    = list(31L, 99L, 3L, 30L, FALSE),
  coll30eq  = list(31L, 99L, 3L, 30L, TRUE)
)

get_t0 <- function(phy, wd) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
             file.path(wd, "cstest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  invisible(suppressWarnings(system2(TNT, args = "cstest.run;", stdout = TRUE, stderr = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  t0
}
run_cfg <- function(phy, t0, cfg, seed) {
  set.seed(seed)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L, nThreads = 1L,
        maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
        xssRounds = 0L, cssRounds = 0L, rssRounds = ROUNDS, wagnerStarts = 1L,
        fuseInterval = 9999L,
        sectorMinSize = cfg[[1]], sectorMaxSize = cfg[[2]], rasStarts = cfg[[3]],
        sectorCollapseTarget = cfg[[4]], sectorAcceptEqual = cfg[[5]]))
  min(as.double(attr(r, "score")))
}

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]]); n <- NTip(phy); tgt <- target[[nm]]
  wd <- file.path(tempdir(), paste0("cs", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  t0 <- get_t0(phy, wd); t0len <- TreeLength(t0, phy)
  cs <- cladeSizes(t0); inband <- sum(cs >= 31 & cs <= 99)
  cat(sprintf("\n==== %s (%dt) | T0=%.0f target=%d | eligible clades in [31,99]: %d  (>30 total: %d) ====\n",
              nm, n, t0len, tgt, inband, sum(cs > 30)))
  if (t0len < tgt - 0.5) cat("  [!] T0 already below target -- mapping/score sanity FAIL; skip\n")
  for (cn in names(cfgs)) {
    sc <- vapply(SEEDS, function(s) run_cfg(phy, t0, cfgs[[cn]], s), double(1))
    best <- min(sc)
    cat(sprintf("  %-10s seeds[%s] -> %s  | best %.0f (%+.0f vs T0, %+.0f vs target)%s\n",
                cn, paste(SEEDS, collapse = ","), paste(format(sc), collapse = " "),
                best, best - t0len, best - tgt, if (best <= tgt) "  <== REACHED" else ""))
  }
}
