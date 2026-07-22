# =====================================================================
# Axis A: is TreeSearch's RAS Wagner consistently DIFFERENT from TNT's?
# Compares the DISTRIBUTION (score + diversity + treespace occupancy) of
# random-addition-sequence Wagner trees (NO branch swapping) between:
#   TSrand   - TreeSearch ts_random_wagner_tree (deterministic first-found
#              tie-break; randomness only from the addition order)
#   TNTdet   - TNT  mult=wagner ras, default rseed]  (deterministic insertion)
#   TNTrand  - TNT  mult=wagner ras, rseed[          (random insertion = random
#              tie-break among equal-best positions: diversity source we lack)
#
# A worse Wagner SCORE is not bad per se (the lead's point: a bad score is
# useful if randomness reaches a basin we'd never find post-TBR).  So we report
# score AND diversity AND cross-set occupancy, and only call the methods
# "consistently different" if TSrand systematically departs from TNTdet.
#
# Env: DS (dataset, default Zanol2014), K (trees per arm, default 60),
#      ARMS (comma list; default all), MDS (1 to save an MDS png).
# =====================================================================
suppressMessages({
  library(TreeSearch, lib.loc = "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem")
  library(TreeTools)
  library(TreeDist)
})
TNT   <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
nm    <- Sys.getenv("DS", "Zanol2014")
K     <- as.integer(Sys.getenv("K", "60"))
seed0 <- 1L

phy  <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
taxa <- names(phy); n <- length(taxa)
cat(sprintf("== Wagner Axis A | %s | n=%d tips | K=%d trees/arm ==\n", nm, n, K))

# --- matrices for the C++ Wagner kernel (mirror AdditionTree.R) ----------
at <- attributes(phy)
contrast <- at$contrast
tipData  <- matrix(unlist(phy, use.names = FALSE), nrow = n, byrow = TRUE)
weight   <- TreeSearch:::.ScaleWeight(at$weight)
levels   <- at$levels

.EdgeToPhylo <- function(edge) {
  tr <- structure(list(edge = edge, tip.label = taxa, Nnode = n - 1L),
                  class = "phylo")
  Renumber(tr)
}

# --- TS random RAS Wagner: K trees, distinct seeds -----------------------
TSrandTrees <- function(k) {
  trees <- vector("list", k)
  for (i in seq_len(k)) {
    set.seed(seed0 + i)
    res <- TreeSearch:::ts_random_wagner_tree(contrast, tipData, weight, levels)
    trees[[i]] <- .EdgeToPhylo(res$edge)
  }
  structure(trees, class = "multiPhylo")
}

# --- TNT RAS Wagner (no swap): K trees ----------------------------------
TNTWagnerTrees <- function(k, randInsert) {
  tag <- if (randInsert) "R" else "D"
  wd  <- file.path(tempdir(), paste0("wag", Sys.getpid(), nm, tag))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c(
    "mxram 1024;",
    "proc data.tnt;",
    sprintf("hold %d;", k + 50L),
    sprintf("rseed %d;", seed0),
    if (randInsert) "rseed[;" else "rseed];",
    sprintf("mult = wagner replic %d keepall;", k),
    "tsave *wag.tre;", "save;", "tsave/;",
    "length;",
    "quit;")
  writeLines(script, file.path(wd, "wagrun.run"))
  old <- setwd(wd)
  out <- suppressWarnings(system2(TNT, args = "wagrun.run;",
                                  stdout = TRUE, stderr = TRUE))
  setwd(old)
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  trees <- ReadTntTree(file.path(wd, "wag.tre"))
  if (!inherits(trees, "multiPhylo")) trees <- structure(list(trees), class = "multiPhylo")
  attr(trees, "tntout") <- out
  trees
}

# --- score every tree with the SAME scorer (TreeLength) ------------------
Scores <- function(trees) vapply(trees, TreeLength, double(1), phy)

# --- run arms ------------------------------------------------------------
armsWanted <- strsplit(Sys.getenv("ARMS", "TSrand,TNTdet,TNTrand"), ",")[[1]]
arms <- list()
if ("TSrand"  %in% armsWanted) arms$TSrand  <- TSrandTrees(K)
if ("TNTdet"  %in% armsWanted) arms$TNTdet  <- TNTWagnerTrees(K, randInsert = FALSE)
if ("TNTrand" %in% armsWanted) arms$TNTrand <- TNTWagnerTrees(K, randInsert = TRUE)

# --- score distributions -------------------------------------------------
cat("\n-- score distribution (TreeLength) --\n")
sc <- lapply(arms, Scores)
for (a in names(arms)) {
  s <- sc[[a]]
  cat(sprintf("  %-8s n=%3d  mean=%.1f sd=%.1f  min=%.0f max=%.0f  distinct.topol=%d\n",
              a, length(s), mean(s), sd(s), min(s), max(s),
              length(unique(arms[[a]]))))
}

# --- within-arm diversity (mean pairwise ClusteringInfoDist) -------------
cat("\n-- within-arm diversity (mean pairwise ClusteringInfoDist, normalized) --\n")
divOf <- function(trees) {
  d <- ClusteringInfoDist(trees, normalize = TRUE)
  c(meanPair = mean(d), medNN = median(apply(as.matrix(d) + diag(Inf, length(trees)), 1, min)))
}
for (a in names(arms)) {
  dv <- divOf(arms[[a]])
  cat(sprintf("  %-8s meanPairwise=%.4f  medianNN=%.4f\n", a, dv["meanPair"], dv["medNN"]))
}

# --- score-distribution tests: TSrand vs TNTdet --------------------------
if (all(c("TSrand", "TNTdet") %in% names(arms))) {
  cat("\n-- TSrand vs TNTdet score-distribution tests --\n")
  w <- suppressWarnings(wilcox.test(sc$TSrand, sc$TNTdet))
  k2 <- suppressWarnings(ks.test(sc$TSrand, sc$TNTdet))
  cat(sprintf("  Mann-Whitney p=%.4g | KS D=%.3f p=%.4g | meanDiff(TSrand-TNTdet)=%+.1f\n",
              w$p.value, k2$statistic, k2$p.value, mean(sc$TSrand) - mean(sc$TNTdet)))
}

# --- cross-set occupancy: TSrand vs TNTdet (within-vs-cross NN) ----------
if (all(c("TSrand", "TNTdet") %in% names(arms))) {
  cat("\n-- treespace occupancy: TSrand vs TNTdet (ClusteringInfoDist) --\n")
  pooled <- structure(c(arms$TSrand, arms$TNTdet), class = "multiPhylo")
  D <- as.matrix(ClusteringInfoDist(pooled, normalize = TRUE))
  k1 <- length(arms$TSrand); idxA <- seq_len(k1); idxB <- k1 + seq_len(length(arms$TNTdet))
  nnIn  <- function(rows, cols) { m <- D[rows, cols, drop = FALSE]; diag(m[, match(rows, cols), drop = FALSE]) <- Inf;
    apply(m, 1, function(r) min(r[is.finite(r)])) }
  # within-set NN (exclude self) and cross-set NN
  withinA <- sapply(idxA, function(i) min(D[i, setdiff(idxA, i)]))
  withinB <- sapply(idxB, function(i) min(D[i, setdiff(idxB, i)]))
  crossAB <- sapply(idxA, function(i) min(D[i, idxB]))
  crossBA <- sapply(idxB, function(i) min(D[i, idxA]))
  cat(sprintf("  within TSrand NN  median=%.4f | within TNTdet NN median=%.4f\n",
              median(withinA), median(withinB)))
  cat(sprintf("  cross TSrand->TNTdet NN median=%.4f | cross TNTdet->TSrand NN median=%.4f\n",
              median(crossAB), median(crossBA)))
  cat("  (cross approx within => same region; cross >> within => different regions)\n")
  if (nzchar(Sys.getenv("MDS"))) {
    pts <- cmdscale(as.dist(D), k = 2)
    png(sprintf("dev/benchmarks/wagner_mds_%s.png", nm), width = 800, height = 800)
    plot(pts, col = c(rep("red", k1), rep("blue", length(arms$TNTdet))), pch = 19,
         main = sprintf("%s Wagner treespace (red=TSrand blue=TNTdet)", nm),
         xlab = "MDS1", ylab = "MDS2")
    dev.off()
    cat(sprintf("  MDS saved: dev/benchmarks/wagner_mds_%s.png\n", nm))
  }
}

# smoke aid: show the TNT length report tail so we can eyeball no-swap behaviour
if (!is.null(arms$TNTdet)) {
  to <- attr(arms$TNTdet, "tntout")
  cat("\n-- TNTdet tail (verify no-swap: lengths vary & are suboptimal) --\n")
  cat(tail(to, 6), sep = "\n"); cat("\n")
}
