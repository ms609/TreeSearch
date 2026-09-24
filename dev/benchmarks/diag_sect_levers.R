# LEVER VERIFICATION (advisor steps 2-3): re-run every sector lever on the now-
# trustworthy harness, reporting BOTH score and candidates_evaluated so we can see
# which flags actually ENGAGE (dCand != 0) vs which are dead wiring (dCand ~ 0).
# Shared start = identical TNT T0 (a TBR-local optimum). Leading hypothesis: our
# default POLISHES the sector (TBR, already stuck); Goloboff's RSS REBUILDS it
# (RAS+TBR, rasStarts>1). If rebuild engages AND drops the score toward TNT, that
# is the missing mechanism. If a lever changes nothing observable, it was never
# tested -- fix the wiring.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006")), "\\s+")[[1]]
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), paste0("levers", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

get_t0 <- function(phy, seed = 1) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", 8), "quit;")
  writeLines(script, file.path(wd, "levtest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "levtest.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, tnt = if (length(s_sect)) s_sect[length(s_sect)] else NA)
}
run <- function(d, tree, rss, ras, ae, mh, ct) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, rssRounds = as.integer(rss),
    rasStarts = as.integer(ras), sectorAcceptEqual = ae,
    sectorMaxHits = as.integer(mh), sectorCollapseTarget = as.integer(ct)))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}
levers <- list(
  base    = list(rss = 8, ras = 1,  ae = FALSE, mh = 1,  ct = 0),
  ras3    = list(rss = 8, ras = 3,  ae = FALSE, mh = 1,  ct = 0),
  ras10   = list(rss = 8, ras = 10, ae = FALSE, mh = 1,  ct = 0),
  ae      = list(rss = 8, ras = 1,  ae = TRUE,  mh = 1,  ct = 0),
  mh20    = list(rss = 8, ras = 1,  ae = FALSE, mh = 20, ct = 0),
  ct10    = list(rss = 8, ras = 1,  ae = FALSE, mh = 1,  ct = 10),
  all     = list(rss = 8, ras = 10, ae = TRUE,  mh = 20, ct = 10))

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  tn <- get_t0(phy)
  start <- TreeLength(tn$t0, phy)
  cat(sprintf("\n==== %s | start(T0)=%.0f  TNT_sect=%.0f  (gap to beat = %+.0f) ====\n",
              nm, start, tn$tnt, tn$tnt - start))
  b <- run(phy, tn$t0, 8, 1, FALSE, 1, 0)
  cat(sprintf("  %-7s score=%.0f cand=%.0f\n", "base", b$score, b$cand))
  for (lv in names(levers)[-1]) {
    p <- levers[[lv]]
    r <- run(phy, tn$t0, p$rss, p$ras, p$ae, p$mh, p$ct)
    cat(sprintf("  %-7s score=%.0f cand=%.0f | dScore=%+.0f dCand=%+.0f %s\n",
                lv, r$score, r$cand, r$score - b$score, r$cand - b$cand,
                ifelse(abs(r$cand - b$cand) < 1, "<-- DEAD (no engage)", "")))
  }
}
