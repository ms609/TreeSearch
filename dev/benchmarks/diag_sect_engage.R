# ENGAGEMENT TEST (advisor steps 1-2): is our sectorial search actually executing,
# and does the rssRounds flag engage? The harness has lied twice today (c.run bug;
# acceptequal==greedy), so the "five nulls" are suspect. From the identical TNT T0
# (a TBR-local optimum), run with sectorial OFF (rssRounds=0) vs ON (rssRounds=8),
# all else fixed, and compare BOTH score and candidates_evaluated.
#   dCand > 0  => sectorial is doing work (executing)
#   dCand ~ 0  => sectorial is NOT running (gated off / dead wiring) -- that's the bug
#   dScore < 0 => sectorial escapes the local optimum (runs AND helps)
#   dScore = 0 => executes but finds nothing (quality/selection/reduction bug)
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Zanol2014")), "\\s+")[[1]]
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), paste0("engage", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

get_t0 <- function(phy, seed = 1) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", 8), "quit;")
  writeLines(script, file.path(wd, "engtest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "engtest.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, tnt = if (length(s_sect)) s_sect[length(s_sect)] else NA)
}
run <- function(d, tree, rss) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, rssRounds = as.integer(rss)))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  tn <- get_t0(phy)
  start <- TreeLength(tn$t0, phy)
  a0 <- run(phy, tn$t0, 0L)
  aK <- run(phy, tn$t0, 8L)
  cat(sprintf("%-11s | start=%.0f TNT_sect=%.0f | rss0: score=%.0f cand=%.0f | rss8: score=%.0f cand=%.0f | dScore=%+.0f dCand=%+.0f\n",
              nm, start, tn$tnt, a0$score, a0$cand, aK$score, aK$cand,
              aK$score - a0$score, aK$cand - a0$cand))
}
