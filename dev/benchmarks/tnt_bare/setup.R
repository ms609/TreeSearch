# Setup: generate Zanol Fitch matrix + RDS, build fixed T0 (TNT mult replic 1, rseed 1).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"), winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }

bare <- "dev/benchmarks/tnt_bare"
nm <- Sys.getenv("DS", "Zanol2014")
phy <- fitch(inapplicable.phyData[[nm]])
WriteTntCharacters(phy, file.path(bare, paste0(nm, ".tnt")))
saveRDS(phy, file.path(bare, paste0(nm, ".phy.rds")))
cat(sprintf("%s: %d tips, %d chars\n", nm, length(phy), attr(phy, "nr")))

# Build a fixed T0 with a single mult replicate (rseed 1). This is our fixture.
wd <- file.path(tempdir(), paste0("t0", Sys.getpid(), nm)); unlink(wd, recursive = TRUE)
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))
writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
             "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;", "quit;"),
           file.path(wd, "buildtee.run"))
old <- setwd(wd)
out <- suppressWarnings(system2(TNT, args = "buildtee.run;", stdout = TRUE, stderr = TRUE))
setwd(old)
file.copy(file.path(wd, "t0.tre"), file.path(bare, paste0(nm, ".t0.tre")), overwrite = TRUE)
t0 <- ReadTntTree(file.path(bare, paste0(nm, ".t0.tre")))
if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
cat(sprintf("T0 (mult replic 1, rseed 1) score = %.0f tips=%d\n", TreeLength(t0, phy), length(t0$tip.label)))
out <- iconv(out, from = "", to = "UTF-8", sub = "")
cat("TNT mult lines:\n"); cat(paste0("  ", grep("score", out, ignore.case = TRUE, value = TRUE)), sep = "\n")
