# Build the SHARED start trees (TNT mult, rseed 1, replic 1) for each gap dataset
# and save as Newick, so the heavy ceiling/validation sweep (local or Hamilton) is
# decoupled from TNT (Hamilton has no TNT).  Also writes the Fitch matrices as
# TNT/phyDat-independent RDS so the sweep needs only TreeSearch + the .tre + .rds.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband2"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN <- c("Zanol2014", "Wortley2006", "Zhu2013", "Giles2015")
outdir <- "dev/benchmarks/t0"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  wd <- file.path(tempdir(), paste0("t0b", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
             file.path(wd, "t0build.run"))
  old <- setwd(wd)
  invisible(suppressWarnings(system2(TNT, args = "t0build.run;", stdout = TRUE, stderr = TRUE)))
  setwd(old)
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  t0len <- TreeLength(t0, phy)
  # Persist tree (Newick) + Fitch phyDat (RDS) so the sweep needs no TNT.
  ape::write.tree(t0, file.path(outdir, paste0(nm, ".tre")))
  saveRDS(phy, file.path(outdir, paste0(nm, ".phy.rds")))
  cat(sprintf("%-12s n=%d  T0 len=%.0f  -> %s.tre + .phy.rds\n",
              nm, NTip(phy), t0len, nm))
}
cat("\nSaved T0 trees + Fitch matrices to", outdir, "\n")
