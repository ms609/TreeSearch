# CHIP FINDING TEST: TNT escapes via sectorial over a DIVERSE SET of equal-optimal
# trees, not single-tree polish.  TreeSearch's rss_search is single-tree-per-replicate
# and MaximizeParsimony(tree=multiPhylo) keeps only tree[[1]] -> it CANNOT operate over
# a set.  Here we test the weaker "independent lanes" route the chip measured (~1/15 reach
# 1261): run our best single-tree sectorial (large-clade [31,99] coll30, 20 picks x 30
# rounds, ratchet off) from EACH of TNT's diverse hold-1000 trees, best-of.  If lanes from
# diverse starts reach the target where the single canonical T0 stalls at the 1267-class
# plateau, start-diversity is confirmed as the lever on our side too.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband2"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m=="-"] <- "?"; MatrixToPhyDat(m) }
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014")), "\\s+")[[1]]
ROUNDS <- as.integer(Sys.getenv("TS_RSSROUNDS", "30"))
SEEDS  <- as.integer(strsplit(Sys.getenv("TS_SEEDS", "1 2 3"), "\\s+")[[1]])
target <- c(Zanol2014 = 1261, Wortley2006 = 480, Zhu2013 = 624, Giles2015 = 670)

# Diverse equal-optimal set from TNT hold-1000 mult (the trees TNT runs sectorial over).
diverse_set <- function(phy, wd) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
               "mult=replic 1;", "tsave *set.tre;", "save;", "tsave/;", "quit;"),
             file.path(wd, "setbuild.run"))
  old <- setwd(wd); on.exit(setwd(old))
  invisible(suppressWarnings(system2(TNT, args = "setbuild.run;", stdout = TRUE, stderr = TRUE)))
  ts <- ReadTntTree(file.path(wd, "set.tre"))
  if (!inherits(ts, "multiPhylo")) ts <- structure(list(ts), class = "multiPhylo")
  ts
}

lane <- function(phy, t, seed) {
  set.seed(seed)
  Sys.setenv(TS_RSS_PICKS = "20")
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t, maxReplicates = 1L, nThreads = 1L,
        maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
        xssRounds = 0L, cssRounds = 0L, rssRounds = ROUNDS, wagnerStarts = 1L,
        fuseInterval = 9999L, sectorMinSize = 31L, sectorMaxSize = 99L,
        rasStarts = 3L, sectorCollapseTarget = 30L, sectorAcceptEqual = FALSE))
  Sys.unsetenv("TS_RSS_PICKS")
  min(as.double(attr(r, "score")))
}

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]]); tgt <- target[[nm]]
  wd <- file.path(tempdir(), paste0("ds", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  ts <- diverse_set(phy, wd)
  lens <- vapply(ts, TreeLength, double(1), phy)
  nset <- length(ts)
  cat(sprintf("\n==== %s | TNT diverse set: %d trees, lengths %.0f-%.0f | target=%d ====\n",
              nm, nset, min(lens), max(lens), tgt))
  allsc <- c()
  for (i in seq_len(nset)) {
    sc <- vapply(SEEDS, function(s) lane(phy, ts[[i]], s), double(1))
    allsc <- c(allsc, sc)
    cat(sprintf("  tree %2d (len %.0f): %s\n", i, lens[i], paste(format(sc), collapse = " ")))
  }
  nhit <- sum(allsc <= tgt + 1e-6)
  cat(sprintf("  >>> %d lanes; best %.0f (target %d); reached target: %d/%d lanes\n",
              length(allsc), min(allsc), tgt, nhit, length(allsc)))
}
