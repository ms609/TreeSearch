suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"), winslash = "/"))
  library(TreeTools)
})
TNT <- "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) suppressWarnings(as.double(gsub(",", "", x)))
nm <- Sys.getenv("DS", "Zanol2014")
phy <- fitch(inapplicable.phyData[[nm]])
wd <- file.path(tempdir(), paste0("seq", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))
writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1;",
             "mult=replic 1;", rep("sectsch=rss;", 8), "quit;"),
           file.path(wd, "seq.run"))
old <- setwd(wd)
out <- suppressWarnings(system2(TNT, args = "seq.run;", stdout = TRUE, stderr = TRUE))
setwd(old)
out <- iconv(out, from = "", to = "UTF-8", sub = "")
cat(sprintf("==== %s | per-pass RSS trace ====\n", nm))
# Show every line that reports a score or replacements, in order
keep <- grep("Best score|RSS|eplac|earrang|ector", out, ignore.case = TRUE, value = TRUE)
cat(paste0("  ", trimws(keep)), sep = "\n")
