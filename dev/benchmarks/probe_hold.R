# RESOLVE the hold-1000-vs-10-trees conflation (user catch).
# `hold 1000` = buffer CAP (1000); the "10 trees" is the emergent count mult deposits.
# Decisive question for a beam design: does TNT's escape need only the ~10 SEED trees, or
# a buffer that keeps GROWING during sectorial?  Isolate by building the SAME 10-tree 1271
# set, then running sectsch under different buffer caps (>=10 so the seed set is intact;
# cap=10 forbids growth, cap=1000 allows it).
suppressMessages({ library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband2"), winslash = "/")); library(TreeTools) })
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m=="-"] <- "?"; MatrixToPhyDat(m) }
nm <- Sys.getenv("TS_DS", "Zanol2014"); phy <- fitch(inapplicable.phyData[[nm]])
wd <- file.path(tempdir(), paste0("ph", Sys.getpid())); unlink(wd, recursive = TRUE)
dir.create(wd, recursive = TRUE, showWarnings = FALSE); WriteTntCharacters(phy, file.path(wd, "data.tnt"))
bestLen <- function(f) { tr <- ReadTntTree(f); if (inherits(tr,"multiPhylo")) min(vapply(tr,TreeLength,double(1),phy)) else TreeLength(tr,phy) }
ntree <- function(f) { tr <- ReadTntTree(f); if (inherits(tr,"multiPhylo")) length(tr) else 1L }

# 1. Build the canonical 10-tree 1271 set (hold 1000 mult), save it.
writeLines(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","mult=replic 1;",
             "tsave *set10.tre;","save;","tsave/;","quit;"), file.path(wd,"build.run"))
old <- setwd(wd); invisible(suppressWarnings(system2(TNT, "build.run;", stdout=TRUE, stderr=TRUE))); setwd(old)
cat(sprintf("==== %s | seed set: %d trees @ len %.0f ====\n", nm, ntree(file.path(wd,"set10.tre")), bestLen(file.path(wd,"set10.tre"))))

# 2. Re-load the SAME 10-tree set, run sectsch under varying buffer cap.
for (cap in c(10L, 12L, 25L, 1000L)) {
  writeLines(c("mxram 1024;","proc data.tnt;","rseed 1;",
               sprintf("hold %d;", cap), "proc set10.tre;",        # load seed set under this cap
               rep("sectsch=rss;", 10), "tsave *out.tre;","save;","tsave/;","quit;"),
             file.path(wd,"sect.run"))
  old <- setwd(wd); invisible(suppressWarnings(system2(TNT, "sect.run;", stdout=TRUE, stderr=TRUE))); setwd(old)
  cat(sprintf("  hold %4d -> land %.0f  (final buffer %d trees)\n",
              cap, bestLen(file.path(wd,"out.tre")), ntree(file.path(wd,"out.tre"))))
}
