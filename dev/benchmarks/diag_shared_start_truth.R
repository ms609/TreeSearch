# DISPOSITIVE shared-start test (advisor): does TNT sectsch reach the target from
# the EXACT T0 our sectorial uses?  Prior "TNT sectsch -> 1261" came from sectsch
# running on TNT's in-memory mult tree, NEVER verified == the t0.tre we fed our
# sectorial -> possible apples-to-oranges.  Here ONE TNT run: mult builds A, saves
# A to t0.tre, runs sectsch FROM A; our sectorial reads the SAME t0.tre.  Both
# share A by construction.  MAPPING CHECK: TreeLength(read t0.tre) must be sane
# (~ the mult score); garbage => ReadTntTree permuted taxa, result invalid.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) suppressWarnings(as.double(gsub(",", "", x)))
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Wortley2006 Zhu2013 Giles2015")), "\\s+")[[1]]
target <- c(Wortley2006 = 479, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  wd <- file.path(tempdir(), paste0("sst", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1000;",
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
               rep("sectsch=rss;", 8), "quit;"), file.path(wd, "ssttest.run"))
  old <- setwd(wd)
  out <- suppressWarnings(system2(TNT, args = "ssttest.run;", stdout = TRUE, stderr = TRUE))
  setwd(old)
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  sect_vals <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                       grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  A_len <- TreeLength(t0, phy)
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L, nThreads = 1L,
        maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L,
        cssRounds = 0L, rssRounds = 8L, rasStarts = 3L, wagnerStarts = 1L, fuseInterval = 9999L))
  ours <- min(as.double(attr(r, "score")))
  tnt_sect <- if (length(sect_vals)) min(sect_vals, na.rm = TRUE) else NA
  cat(sprintf("\n==== %s | target=%d ====\n", nm, target[[nm]]))
  cat(sprintf("  TNT mult T0 (our ruler) A_len = %.0f   [mapping sane? expect a real MP-ish score]\n", A_len))
  cat(sprintf("  TNT sectsch FROM A      -> %s   (escape %+.0f vs A)\n",
              format(tnt_sect), if (is.finite(tnt_sect)) tnt_sect - A_len else NA))
  cat(sprintf("  OUR sectorial FROM A    -> %.0f   (escape %+.0f vs A)\n", ours, ours - A_len))
  cat(sprintf("  VERDICT: %s\n",
              if (is.finite(tnt_sect) && tnt_sect < A_len - 0.5)
                "TNT sectsch ESCAPES shared A -> real sectorial gap, trace mechanism"
              else "TNT sectsch does NOT escape shared A -> 1261 was a different basin; hunt dissolves"))
}
