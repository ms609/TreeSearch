# Fidelity oracle (advisor): CSS runs sector-restricted TBR against the FULL
# dataset (ts_sector.cpp css_search) -> EXACT scoring, no HTU approximation.
# If CSS made to search hard STILL can't match TNT from an identical T0, the HTU
# fidelity hypothesis is exonerated and the gap is STRUCTURAL (sector shape /
# selection). If exact-scoring CSS closes it where HTU-based RSS could not,
# fidelity is the culprit.
#
# Shared-start design (as bench_ras_verify.R). Arms from TNT's own T0:
#   rss   = HTU-based RSS, K passes        (the gap baseline)
#   css2  = exact CSS, K rounds, 2 big sectors
#   css3  = exact CSS, K rounds, 3 sectors
# gap = TS_sect - TNT_sect; lower = closer. Also reports Mcand (did CSS search?).
#
# Env: TS_LIB (default .agent-ratchet), TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS, OUT_CSV.
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
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/sect_css_oracle.csv")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "cssoracle"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

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
run_rss <- function(d, tree, rounds) {
  set.seed(1); nt <- length(d)
  smin <- as.integer(round(nt * 0.35)); smax <- as.integer(round(nt * 0.65))
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, sectorMinSize = smin, sectorMaxSize = smax,
    rssRounds = as.integer(rounds)))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}
run_css <- function(d, tree, rounds, part) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, rssRounds = 0L,
    cssRounds = as.integer(rounds), cssPartitions = as.integer(part),
    wagnerStarts = 1L, fuseInterval = 9999L))
  list(score = as.double(attr(r, "score")), cand = as.double(attr(r, "candidates_evaluated")))
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s s%d: no T0\n", nm, sd)); next }
    rss  <- run_rss(phy, tn$t0, K)
    css2 <- run_css(phy, tn$t0, K, 2L)
    css3 <- run_css(phy, tn$t0, K, 3L)
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd, tnt = tn$s_sect,
      rss = rss$score, css2 = css2$score, css3 = css3$score,
      g_rss = rss$score - tn$s_sect, g_css2 = css2$score - tn$s_sect,
      g_css3 = css3$score - tn$s_sect,
      Mc_rss = round(rss$cand / 1e6, 1), Mc_css2 = round(css2$cand / 1e6, 1),
      stringsAsFactors = FALSE)
    cat(sprintf("%-11s s%d | TNT=%.0f | rss=%.0f css2=%.0f css3=%.0f | gaps %+.0f/%+.0f/%+.0f | Mc rss=%.1f css2=%.1f\n",
                nm, sd, tn$s_sect, rss$score, css2$score, css3$score,
                rss$score - tn$s_sect, css2$score - tn$s_sect, css3$score - tn$s_sect,
                rss$cand / 1e6, css2$cand / 1e6))
  }
}
S <- do.call(rbind, rows)
cat("\n== medians (gap = TS_sect - TNT_sect from identical T0; lower = closer) ==\n")
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], TNT = median(d$tnt),
  g_rss = median(d$g_rss), g_css2 = median(d$g_css2), g_css3 = median(d$g_css3),
  Mc_rss = median(d$Mc_rss), Mc_css2 = median(d$Mc_css2))))
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
