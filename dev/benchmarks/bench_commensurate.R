# COMMENSURABILITY CHECK (advisor): is TNT's stdout sectorial score on the same
# scale as our TreeLength? The shared-start "gap" compared our TreeLength-scored
# tree to TNT's STDOUT number for a tree we never re-scored. Here we save TNT's
# actual best tree and score it ourselves.
#   TreeLength(TNT_best) ~ TreeLength(T0) (e.g. both 1275) while TNT stdout says
#     T0=1275 sect=1262  => SCORING OFFSET, same tree, NO search gap.
#   TreeLength(TNT_best) genuinely < TreeLength(T0) => real improvement, real gap.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Wortley2006 Zhu2013 Giles2015")), "\\s+")[[1]]
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
tl_any <- function(tr, phy) {
  if (is.null(tr)) return(NA_real_)
  if (inherits(tr, "multiPhylo")) min(vapply(tr, TreeLength, double(1), phy)) else TreeLength(tr, phy)
}
# Unique per-process temp dir: a stale/locked data.tnt from a dead TNT orphan in
# a shared dir makes the next run fail silently (NA). Script stem MUST be purely
# alphabetic and not a TNT command -- TNT parses the filename as its command line
# (see dev/expertise/tnt.md). "c.run" => command `c` (ccode) => "Must read data
# before changing character settings"; "cmnstest.run" is safe.
wd <- file.path(tempdir(), paste0("commens", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", "rseed 1;", "taxname=;",
              "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", 8), "tsave *best.tre;", "save;", "tsave/;", "quit;")
  writeLines(script, file.path(wd, "cmnstest.run"))
  old <- setwd(wd)
  out <- suppressWarnings(system2(TNT, args = "cmnstest.run;", stdout = TRUE, stderr = TRUE))
  setwd(old)
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  if (!any(grepl("Best score", out))) {
    cat(sprintf("==== %s: TNT produced NO score; raw output ====\n", nm))
    cat(head(out, 30), sep = "\n"); cat("\n")
  }
  s_tbr <- num(sub(".*Best score \\(TBR\\):\\s*([0-9.]+).*", "\\1",
                   grep("Best score \\(TBR\\):", out, value = TRUE)[1]))
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  s_sect <- if (length(s_sect)) s_sect[length(s_sect)] else NA
  t0   <- tryCatch(ReadTntTree(file.path(wd, "t0.tre")),   error = function(e) NULL)
  best <- tryCatch(ReadTntTree(file.path(wd, "best.tre")), error = function(e) NULL)
  cat(sprintf("%-11s | TNT stdout: T0=%.0f sect=%.0f | TreeLength: T0=%.0f best=%.0f\n",
              nm, s_tbr, s_sect, tl_any(t0, phy), tl_any(best, phy)))
}
