suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"), winslash = "/"))
  library(TreeTools)
})
TNT <- "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) suppressWarnings(as.double(gsub(",", "", x)))
dsN <- c("Wortley2006", "Zanol2014", "Zhu2013", "Giles2015")
target <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
cfgs <- list(
  default  = character(0),
  noglobal = "sectsch: noglobal;",
  equals   = "sectsch: equals;",
  global1  = "sectsch: global 1;",
  recurse2 = "sectsch: recurse 2;"
)
rx_best <- "Sectorial search \\(RSS\\), best score:"
rx_tbr  <- "Best score \\(TBR\\):"
run_cfg <- function(wd, setlines) {
  is_rec <- any(grepl("recurse", setlines))
  pre  <- if (is_rec) setlines else character(0)
  post <- if (is_rec) character(0) else setlines
  writeLines(c("mxram 1024;", pre, "proc data.tnt;", "rseed 1;", "hold 1;",
               "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
               post, rep("sectsch=rss;", 8), "quit;"),
             file.path(wd, "optest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "optest.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  vl <- grep(rx_best, out, value = TRUE)
  v <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1", vl))
  tl <- grep(rx_tbr, out, value = TRUE)
  t0 <- num(sub(".*\\(TBR\\):\\s*([0-9.]+).*", "\\1", tl[1]))
  list(t0 = t0, best = if (length(v)) min(v, na.rm = TRUE) else NA_real_)
}
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  wd <- file.path(tempdir(), paste0("ng", Sys.getpid(), nm))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  cat(sprintf("\n==== %s | target(sectsch)=%d ====\n", nm, target[[nm]]))
  for (cn in names(cfgs)) {
    r <- run_cfg(wd, cfgs[[cn]])
    d <- if (is.finite(r$best) && is.finite(r$t0)) r$best - r$t0 else NA_real_
    cat(sprintf("  %-9s T0=%s sectsch_best=%s (%+.0f vs T0)\n",
                cn, format(r$t0), format(r$best), d))
  }
}
