# WORTLEY 480->479 ESCAPE PROBE — the simplest reproducible instance of the gap.
# TS stalls at 480 (15/15 runs, even intensive+fuse, 60s); TNT xmult reaches 479.
# Characterise the barrier on the SMALLEST case:
#   (a) commensurate: does TreeLength score TNT's 479 tree as 479?
#   (b) are 479 and TS's 480 BOTH TS-TBR local optima?  (if 479 drops below under
#       our TBR, 479 isn't optimal; if both hold, they are separate basins)
#   (c) how far apart are the two optima? (RF split distance)
#   (d) can TS's TBR reach 479 if STARTED adjacent to it? (perturb T479 by 1 NNI,
#       re-search: returns to 479 => findable basin; falls to 480 => leaks away)
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
phy <- fitch(inapplicable.phyData[["Wortley2006"]])

# RF (symmetric split difference); RobinsonFoulds not exported in this TreeTools.
rf <- function(a, b) {
  b <- KeepTip(b, a$tip.label)
  cl <- function(t) {
    pp <- ape::prop.part(t); labs <- attr(pp, "labels")
    s <- vapply(pp, function(ix) {
      x <- labs[ix]; if (length(x) > length(labs) / 2) x <- setdiff(labs, x)
      paste(sort(x), collapse = ",")
    }, character(1))
    s[vapply(strsplit(s, ","), length, 1L) >= 2L]
  }
  sa <- cl(a); sb <- cl(b)
  length(setdiff(sa, sb)) + length(setdiff(sb, sa))
}

# Pure-TBR polish from a given tree (ratchet/drift/sectorial OFF).
polish <- function(tree) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = tree, maxReplicates = 1L,
        nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
        ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, rssRounds = 0L,
        cssRounds = 0L, wagnerStarts = 1L, fuseInterval = 9999L))
  min(as.double(attr(r, "score")))
}

# 1. TNT -> 479 tree
wd <- file.path(tempdir(), paste0("esc", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))
writeLines(c("mxram 1024;", "proc data.tnt;", "hold 10000;", "rseed 1;",
             "xmult=hits 10 replic 50;", "best;", "tsave *t479.tre;", "save;",
             "tsave/;", "quit;"), file.path(wd, "esctest.run"))
old <- setwd(wd)
out <- suppressWarnings(system2(TNT, args = "esctest.run;", stdout = TRUE, stderr = TRUE))
setwd(old)
out <- iconv(out, from = "", to = "UTF-8", sub = "")
tnt_score <- num(sub(".*Best score:\\s*([0-9.]+).*", "\\1", grep("Best score:", out, value = TRUE)[1]))
T479 <- ReadTntTree(file.path(wd, "t479.tre")); if (inherits(T479, "multiPhylo")) T479 <- T479[[1]]
len479 <- TreeLength(T479, phy)
cat(sprintf("(a) TNT Best score=%.0f | TreeLength(T479)=%.0f  [commensurate: %s]\n",
            tnt_score, len479, if (isTRUE(tnt_score == len479)) "YES" else "NO!"))

# 2. TS from scratch -> best (intensive, 3 seeds)
best_ts <- NULL; best_len <- Inf
for (s in 1:3) {
  set.seed(s)
  r <- suppressWarnings(MaximizeParsimony(phy, strategy = "intensive",
        maxReplicates = 9999L, maxSeconds = 30, nThreads = 1L, verbosity = 0L))
  l <- min(as.double(attr(r, "score")))
  tr <- if (inherits(r, "multiPhylo")) r[[1]] else r
  if (l < best_len) { best_len <- l; best_ts <- tr }
}
cat(sprintf("    TS-from-scratch best=%.0f\n", best_len))

# 3/(b) local-optimum check
p479 <- polish(T479); p480 <- polish(best_ts)
cat(sprintf("(b) TBR-polish T479 -> %.0f (479 %s TBR-optimal) | TBR-polish TS-best -> %.0f (%.0f %s TBR-optimal)\n",
            p479, if (p479 >= len479) "IS" else "NOT", p480, best_len,
            if (p480 >= best_len) "IS" else "NOT"))

# 4/(c) basin distance
cat(sprintf("(c) RF(T479, TS-best) = %d splits (max %d)\n", rf(T479, best_ts), NTip(phy) - 3L))

# 5/(d) start TS adjacent to T479: does it find 479?
set.seed(7)
adj <- TBRMoves <- NULL
adj <- tryCatch(TreeTools::Postorder(ape::rNNI(T479, moves = 1L)), error = function(e) T479)
r_adj <- suppressWarnings(MaximizeParsimony(phy, tree = adj, strategy = "intensive",
          maxReplicates = 50L, maxSeconds = 30, nThreads = 1L, verbosity = 0L))
cat(sprintf("(d) TS from (T479 + 1 NNI) -> %.0f  (%s recover 479)\n",
            min(as.double(attr(r_adj, "score"))),
            if (min(as.double(attr(r_adj, "score"))) <= len479) "DID" else "did NOT"))
