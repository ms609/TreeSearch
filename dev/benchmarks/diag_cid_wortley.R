# CID measurement (owed): how far is TS's Wortley optimum from TNT's 479 tree,
# by TreeDist::ClusteringInfoDist(normalize=TRUE) -- robust, unlike the RF=54 I
# reported (RF inflates when one tip moves far).  Low CID => basins are actually
# close (RF was misleading); high CID => genuine whole-tree difference.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
  library(TreeDist)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy <- fitch(inapplicable.phyData[["Wortley2006"]])

# TNT 479 tree (full xmult)
wd <- file.path(tempdir(), paste0("cid", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))
writeLines(c("mxram 1024;", "proc data.tnt;", "hold 10000;", "rseed 1;",
             "xmult=hits 10 replic 50;", "best;", "tsave *t479.tre;", "save;",
             "tsave/;", "quit;"), file.path(wd, "cidtest.run"))
old <- setwd(wd)
invisible(suppressWarnings(system2(TNT, args = "cidtest.run;", stdout = TRUE, stderr = TRUE)))
setwd(old)
T479 <- ReadTntTree(file.path(wd, "t479.tre")); if (inherits(T479, "multiPhylo")) T479 <- T479[[1]]

# TS optimum (intensive, best of 3 short seeds)
best <- NULL; bestlen <- Inf
for (s in 1:3) {
  set.seed(s)
  r <- suppressWarnings(MaximizeParsimony(phy, strategy = "intensive",
        maxReplicates = 9999L, maxSeconds = 20, nThreads = 1L, verbosity = 0L))
  l <- min(as.double(attr(r, "score")))
  tr <- if (inherits(r, "multiPhylo")) r[[1]] else r
  if (l < bestlen) { bestlen <- l; best <- tr }
}

T479 <- KeepTip(T479, best$tip.label)
cid  <- TreeDist::ClusteringInfoDist(T479, best, normalize = TRUE)
rf   <- TreeDist::RobinsonFoulds(T479, best, normalize = TRUE)
cat(sprintf("TNT T479 len=%.0f | TS-best len=%.0f\n", TreeLength(T479, phy), bestlen))
cat(sprintf("ClusteringInfoDist(normalize=TRUE) = %.3f   [0=identical, 1=maximally different]\n", cid))
cat(sprintf("RobinsonFoulds(normalize=TRUE)      = %.3f   (for contrast)\n", rf))
