# DEFINE THE TARGET (advisor's blocking point): from the CANONICAL hold-1000 T0,
# where does TNT's ratchet/drift-OFF RSS sectorial land?  "Did TreeSearch match" is
# undefined until this number exists for the SAME start tree the TS arms will use.
#
# Recipe = ratchet-off sectorial: plain `mult=replic 1` (one RAS+TBR start, NO
# ratchet/drift/fuse) -> T0; then repeated `sectsch=rss;` (each = one full RSS round,
# escape-doc-verified).  hold 1000 selects the canonical 1271 basin (subchip: hold 1
# -> 1275, hold 1000 -> 1271 for Zanol).  Saves each T0 to dev/benchmarks/t0/<nm>.tre
# (Newick) + the Fitch phyDat (.phy.rds) so the TreeSearch arms share the identical start.
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
NROUND <- as.integer(Sys.getenv("TS_NROUND", "10"))
outdir <- "dev/benchmarks/t0"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bestLen <- function(tr, phy) {
  if (inherits(tr, "multiPhylo")) min(vapply(tr, TreeLength, double(1), phy)) else TreeLength(tr, phy)
}

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]]); n <- NTip(phy)
  wd <- file.path(tempdir(), paste0("tgt", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
               "mult=replic 1;",
               "tsave *t0.tre;", "save;", "tsave/;",            # canonical T0 (1271 basin)
               rep("sectsch=rss;", NROUND),                     # ratchet-off RSS rounds
               "tsave *final.tre;", "save;", "tsave/;",         # post-sectorial
               "quit;"),
             file.path(wd, "deftgt.run"))
  old <- setwd(wd)
  out <- suppressWarnings(system2(TNT, args = "deftgt.run;", stdout = TRUE, stderr = TRUE))
  setwd(old)
  out <- iconv(out, from = "", to = "UTF-8", sub = "")

  t0  <- ReadTntTree(file.path(wd, "t0.tre"))
  fin <- ReadTntTree(file.path(wd, "final.tre"))
  t0len  <- bestLen(t0, phy)
  finlen <- bestLen(fin, phy)
  # Persist canonical T0 + matrix for the shared-start TreeSearch arms
  t0one <- if (inherits(t0, "multiPhylo")) t0[[1]] else t0
  ape::write.tree(t0one, file.path(outdir, paste0(nm, ".tre")))
  saveRDS(phy, file.path(outdir, paste0(nm, ".phy.rds")))

  cat(sprintf("%-12s n=%3d | T0=%.0f  -> TNT ratchet-off sectorial(%d rounds) = %.0f  (escape %+.0f)\n",
              nm, n, t0len, NROUND, finlen, finlen - t0len))
}
cat(sprintf("\nSaved canonical T0 trees + matrices to %s/\n", outdir))
