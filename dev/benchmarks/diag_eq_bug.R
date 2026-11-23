# IS sectorAcceptEqual=TRUE BUGGY?  Decisive no-op check.
# Static read (ts_sector.cpp): the equal branch (1181-1187) KEEPS the equal-score
# topology (no revert), so accept_equal is structurally LIVE -- BUT unlike the strict
# branch (1147-1160) it does NOT recompute subtree_size/eligible, so selection goes
# STALE across a plateau walk (latent defect, flagged separately).
# Empirical no-op test: from a near-optimal TNT T0, run rss (rasStarts=3, large
# sectors) with accept_equal F vs T and compare the OUTPUT TOPOLOGY (CID):
#   EQUAL-keep>0 (TS_SECT_DEBUG) AND CID(treeF,treeT)>0 -> LIVE (plateau-walks; equal
#     final score just means no downhill exit) -> NOT a no-op bug.
#   EQUAL-keep==0 -> equal moves never PROPOSED (search_sector tie path dead) -> BUG.
#   EQUAL-keep>0 but CID==0 -> kept-but-identical reinsert -> BUG.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband2"),
            winslash = "/"))
  library(TreeTools)
  library(TreeDist)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
nm <- Sys.getenv("TS_DS", "Zanol2014")
SMIN <- as.integer(Sys.getenv("TS_SMIN", "31"))
SMAX <- as.integer(Sys.getenv("TS_SMAX", "99"))
MAXHITS <- as.integer(Sys.getenv("TS_MAXHITS", "1"))
ROUNDS <- as.integer(Sys.getenv("TS_ROUNDS", "5"))
phy <- fitch(inapplicable.phyData[[nm]])

wd <- file.path(tempdir(), paste0("eqb", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))
writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
             "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
           file.path(wd, "eqbtest.run"))
old <- setwd(wd); invisible(suppressWarnings(system2(TNT, args = "eqbtest.run;", stdout = TRUE, stderr = TRUE))); setwd(old)
t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
t0len <- TreeLength(t0, phy)
cat(sprintf("==== %s | T0 len=%.0f (rasStarts=3, sectors[%d,%d], rssRounds=5) ====\n", nm, t0len, SMIN, SMAX))

Sys.setenv(TS_SECT_DEBUG = "1")   # streams STRICT / EQUAL-keep / WORSE-revert to stderr
run1 <- function(eq) {
  set.seed(1)
  suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L, nThreads = 1L,
    maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
    xssRounds = 0L, cssRounds = 0L, rssRounds = ROUNDS, wagnerStarts = 1L, fuseInterval = 9999L,
    sectorMinSize = SMIN, sectorMaxSize = SMAX, rasStarts = 3L, sectorMaxHits = MAXHITS,
    sectorCollapseTarget = 0L, sectorAcceptEqual = eq))
}
cat("\n--- run F (accept_equal=FALSE) ---\n"); treeF <- run1(FALSE)
cat("\n--- run T (accept_equal=TRUE) ---\n");  treeT <- run1(TRUE)
sF <- min(as.double(attr(treeF, "score"))); sT <- min(as.double(attr(treeT, "score")))
if (inherits(treeF, "multiPhylo")) treeF <- treeF[[1]]
if (inherits(treeT, "multiPhylo")) treeT <- treeT[[1]]
cid <- function(a, b) as.double(ClusteringInfoDist(a, b, normalize = TRUE))
cat(sprintf("\n  accept_equal=F: score=%.0f  CID(T0,F)=%.3f\n", sF, cid(t0, treeF)))
cat(sprintf("  accept_equal=T: score=%.0f  CID(T0,T)=%.3f\n", sT, cid(t0, treeT)))
dFT <- cid(treeF, treeT)
cat(sprintf("  CID(treeF, treeT) = %.3f  %s\n", dFT,
            if (dFT > 1e-6) "<- topology DIFFERS (accept_equal LIVE)" else
              "<- IDENTICAL (accept_equal had NO topological effect)"))
