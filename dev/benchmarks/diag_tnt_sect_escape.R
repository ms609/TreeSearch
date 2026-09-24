# FOUNDATIONAL CHECK for the D1 hunt: does TNT's sectorial (sectsch=rss) actually
# ESCAPE its own mult T0?  The audit's whole premise is "TNT RSS improves T0 by
# +3..+11; ours improves 0".  If TNT sectsch does NOT beat its mult T0, the
# sectorial-escape story is a misattribution and D1 is moot.  TNT-only.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) suppressWarnings(as.double(gsub(",", "", x)))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Zanol2014")), "\\s+")[[1]]

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  wd <- file.path(tempdir(), paste0("se", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  writeLines(c("mxram 1024;", "proc data.tnt;", "hold 1000;", "rseed 1;",
               "mult=replic 1;", "best;",            # T0 (mult) best score
               rep("sectsch=rss;", 8), "best;",       # post-sectorial best score
               "quit;"), file.path(wd, "setest.run"))
  old <- setwd(wd)
  out <- suppressWarnings(system2(TNT, args = "setest.run;", stdout = TRUE, stderr = TRUE))
  setwd(old)
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  best_lines <- grep("Best score:", out, value = TRUE)
  best_vals  <- num(sub(".*Best score:\\s*([0-9.]+).*", "\\1", best_lines))
  sect_lines <- grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)
  sect_vals  <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1", sect_lines))
  t0    <- if (length(best_vals)) best_vals[1] else NA
  final <- if (length(best_vals)) min(best_vals, na.rm = TRUE) else NA
  cat(sprintf("%-11s | mult T0=%s | final=%s | escape=%s | sectsch progression: %s\n",
              nm, format(t0), format(final),
              if (is.finite(t0) && is.finite(final)) sprintf("%+.0f", final - t0) else "NA",
              paste(format(sect_vals), collapse=" ")))
}
