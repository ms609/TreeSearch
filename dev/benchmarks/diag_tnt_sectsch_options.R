# TNT sectsch OPTION TRACE: from the shared T0 (mult, rseed 1), run sectsch=rss
# with escape-relevant knobs toggled, to isolate WHICH drives TNT's -10 escape.
#   noglobal  -> if escape dies, the GLOBAL-TBR cadence is the mechanism
#   equals    -> if escape grows/changes, LATERAL acceptance matters
#   global 1  -> max global-TBR cadence
# TNT-only, deterministic T0 across runs (same rseed/mult).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) suppressWarnings(as.double(gsub(",", "", x)))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014")), "\\s+")[[1]]

# each config = the "sectsch: <set>;" line(s) before the runs ("" = defaults)
cfgs <- list(
  default   = character(0),
  noglobal  = "sectsch: noglobal;",
  equals    = "sectsch: equals;",
  global1   = "sectsch: global 1;",
  eq_global1= c("sectsch: equals;", "sectsch: global 1;")
)
run_cfg <- function(phy, wd, setlines) {
  writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
               "mult=replic 1;", setlines, rep("sectsch=rss;", 8), "quit;"),
             file.path(wd, "optest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "optest.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  v <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
               grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  if (length(v)) min(v, na.rm = TRUE) else NA_real_
}
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  wd <- file.path(tempdir(), paste0("opt", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  cat(sprintf("\n==== %s (T0 ~ TNT mult) | sectsch option trace ====\n", nm))
  for (cn in names(cfgs)) {
    sc <- run_cfg(phy, wd, cfgs[[cn]])
    cat(sprintf("  %-11s sectsch best = %s\n", cn, format(sc)))
  }
}
