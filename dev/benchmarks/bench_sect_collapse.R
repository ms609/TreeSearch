# Does COLLAPSING sub-clades into composite terminals (Goloboff 1999's reduced
# dataset) close the shared-start sectorial gap? sectorCollapseTarget>0 prunes a
# big selected clade to ~target composite first-pass terminals, so the sector
# search rearranges the coarse skeleton of major sub-clades instead of shuffling
# tips within a contiguous clade. Pairs with rasStarts (RAS-rebuild the skeleton
# = TNT's sectsch). Shared-start design; gap = TS_sect - TNT_sect, lower = closer.
#
# Env: TS_LIB (default .agent-ratchet), TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS,
#      TS_COLLAPSE (target terminals, default 10), OUT_CSV.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
K   <- as.integer(Sys.getenv("TS_KPASS", "8"))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
         "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
CT  <- as.integer(Sys.getenv("TS_COLLAPSE", "10"))
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/sect_collapse.csv")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "sectcoll"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

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
run_ts <- function(d, tree, rss, ras, collapse) {
  set.seed(1); nt <- length(d)
  smin <- as.integer(round(nt * 0.35)); smax <- as.integer(round(nt * 0.65))
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, sectorMinSize = smin, sectorMaxSize = smax,
    rssRounds = as.integer(rss), rasStarts = as.integer(ras),
    sectorCollapseTarget = as.integer(collapse)))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s s%d: no T0\n", nm, sd)); next }
    base    <- run_ts(phy, tn$t0, K, 1L, 0L)    # no collapse (gap baseline)
    coll    <- run_ts(phy, tn$t0, K, 1L, CT)    # collapse + TBR
    collras <- run_ts(phy, tn$t0, K, 3L, CT)    # collapse + RAS-rebuild (= TNT)
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd, tnt = tn$s_sect,
      base = base$score, coll = coll$score, collras = collras$score,
      g_base = base$score - tn$s_sect, g_coll = coll$score - tn$s_sect,
      g_collras = collras$score - tn$s_sect,
      Mc_base = round(base$cand / 1e6, 1), Mc_collras = round(collras$cand / 1e6, 1),
      stringsAsFactors = FALSE)
    cat(sprintf("%-11s s%d | TNT=%.0f | base=%.0f coll=%.0f collras=%.0f | gaps %+.0f/%+.0f/%+.0f | Mc %.1f->%.1f\n",
                nm, sd, tn$s_sect, base$score, coll$score, collras$score,
                base$score - tn$s_sect, coll$score - tn$s_sect, collras$score - tn$s_sect,
                base$cand / 1e6, collras$cand / 1e6))
  }
}
S <- do.call(rbind, rows)
cat(sprintf("\n== medians (gap = TS_sect - TNT_sect from identical T0; collapse target=%d) ==\n", CT))
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], TNT = median(d$tnt),
  g_base = median(d$g_base), g_coll = median(d$g_coll), g_collras = median(d$g_collras))))
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
