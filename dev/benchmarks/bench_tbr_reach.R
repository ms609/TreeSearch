# Is TNT's sectorial improvement even TBR-reachable from T0? Five sector mechanisms
# are null; exact-CSS showed T0 is TBR-optimal for us at max_hits=1. This isolates
# whether a THOROUGH global TBR (hold many equal trees -> traverse plateaus) from
# the identical T0 reaches TNT's sectorial score. ratchet/drift/sectors all OFF.
#   tbr50 reaches TNT  => improvement was a global plateau our max_hits=1 TBR missed
#   tbr50 ~ start      => not TBR-reachable; needs rebuild (and our RAS is failing)
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
K   <- as.integer(Sys.getenv("TS_KPASS", "8"))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Zanol2014")), "\\s+")[[1]]
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/tbr_reach.csv")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "tbrreach"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

run_tnt <- function(phy, seed, kpass) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", kpass), "quit;")
  writeLines(script, file.path(wd, "ss.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "ss.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- tryCatch(ReadTntTree(file.path(wd, "t0.tre")), error = function(e) NULL)
  if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, s_sect = if (length(s_sect)) s_sect[length(s_sect)] else NA)
}
run_tbr <- function(d, tree, maxhits) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, rssRounds = 0L,
    cssRounds = 0L, wagnerStarts = 1L, fuseInterval = 9999L,
    tbrMaxHits = as.integer(maxhits)))
  as.double(attr(r, "score"))
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s s%d: no T0\n", nm, sd)); next }
    start <- TreeLength(tn$t0, phy)
    t1  <- run_tbr(phy, tn$t0, 1L)
    t50 <- run_tbr(phy, tn$t0, 50L)
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd, start = start,
      tnt = tn$s_sect, tbr1 = t1, tbr50 = t50,
      g_tbr1 = t1 - tn$s_sect, g_tbr50 = t50 - tn$s_sect, stringsAsFactors = FALSE)
    cat(sprintf("%-11s s%d | start=%.0f TNT=%.0f | tbr1=%.0f tbr50=%.0f | g_tbr50=%+.0f\n",
                nm, sd, start, tn$s_sect, t1, t50, t50 - tn$s_sect))
  }
}
S <- do.call(rbind, rows)
cat("\n== medians (gap = TS_TBR - TNT_sect from identical T0) ==\n")
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], start = median(d$start), TNT = median(d$tnt),
  tbr1 = median(d$tbr1), tbr50 = median(d$tbr50),
  g_tbr50 = median(d$g_tbr50))))
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
