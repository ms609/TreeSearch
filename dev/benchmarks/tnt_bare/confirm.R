# Cross-dataset confirmation of the single-strict-plateaus / set-strict-escapes pattern.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"), winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) suppressWarnings(as.double(gsub(",", "", x)))
target <- c(Zanol2014 = 1261, Wortley2006 = 479, Giles2015 = 670, Zhu2013 = 624)
SEEDS <- 1:4; RDS <- 30

for (nm in strsplit(Sys.getenv("DSETS", "Zanol2014 Wortley2006 Giles2015"), "\\s+")[[1]]) {
  phy <- fitch(inapplicable.phyData[[nm]])
  wd <- file.path(tempdir(), paste0("cf", Sys.getpid(), nm)); unlink(wd, recursive = TRUE)
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  run_tnt <- function(lines) { writeLines(lines, file.path(wd, "runme.run"))
    old <- setwd(wd); o <- suppressWarnings(system2(TNT, "runme.run;", stdout = TRUE, stderr = TRUE))
    setwd(old); iconv(o, from = "", to = "UTF-8", sub = "") }
  best <- function(lines) min(num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
            grep("Sectorial search \\(RSS\\), best score:", run_tnt(lines), value = TRUE))))
  # Build T0 SET (hold 1000) and a single T0 file
  run_tnt(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","mult=replic 1;",
            "tsave *set.tre;","save;","tsave/;","quit;"))
  L <- readLines(file.path(wd, "set.tre"))
  writeLines(c(L[1], paste0(sub("[*]$","",L[2]),";"), "proc-;"), file.path(wd, "tee.tre"))
  start_n <- length(grep("[*]", L)) + 1L
  ss   <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
            "hold 1000;","proc tee.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",RDS),"quit;")))
  se   <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
            "hold 1000;","proc tee.tre;","sectsch: equals;", rep("sectsch=rss;",RDS),"quit;")))
  set  <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
            "hold 1000;","proc set.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",RDS),"quit;")))
  setd <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
            "hold 1000;","proc set.tre;", rep("sectsch=rss;",RDS),"quit;")))  # DEFAULT (noequals) on set
  f <- function(v) sprintf("med=%g [%g-%g]", median(v), min(v), max(v))
  cat(sprintf("\n==== %s (target %d, %d-tree start set) ====\n", nm, target[[nm]], start_n))
  cat(sprintf("  SINGLE-T0 strict : %s\n", f(ss)))
  cat(sprintf("  SINGLE-T0 equals : %s\n", f(se)))
  cat(sprintf("  SET strict       : %s\n", f(set)))
  cat(sprintf("  SET default(TNT) : %s\n", f(setd)))
}
