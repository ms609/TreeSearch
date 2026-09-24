# Diagnostic: characterise the FIRST score-improving TNT sectorial move, to
# decide WHICH fix our sectorial needs (advisor's 4-way discrimination):
#   - non-clade band            -> full multi-stub reduced-dataset rewrite
#   - clade OUTSIDE size band    -> just widen sector selection (trivial)
#   - clade in-band, ATTACHMENT-only change -> b=1 floating-HTU (deferred piece)
#   - clade in-band, INTERNAL change         -> fix RAS diversity/acceptance
#
# Method (advisor): isolate ONE operation. Run TNT `mult` -> T0, then K single
# `sectsch=rss` passes, saving the tree after EACH pass. Diff the FIRST pass that
# drops the score (A = pre, B = post). K=8 cumulative would look non-clade.
#
# Env: TS_LIB, TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 3")), "\\s+")[[1]])
K   <- as.integer(Sys.getenv("TS_KPASS", "8"))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Wortley2006 Zhu2013")),
                "\\s+")[[1]]
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "sectshape"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

# Canonical key per split: sorted tip labels of the side not containing tip[1]
# (size tie broken lexicographically). RF == 0  <=>  equal key sets.
split_keys <- function(tree) {
  sp <- as.Splits(tree)
  m <- as.logical(sp)                       # n_split x n_tip (TRUE = in split)
  if (is.null(dim(m))) m <- matrix(m, nrow = 1)
  labs <- TipLabels(tree)
  apply(m, 1, function(r) {
    a <- sort(labs[r]); b <- sort(labs[!r])
    side <- if (length(a) != length(b)) {
      if (length(a) < length(b)) a else b
    } else if (paste(a, collapse = ",") < paste(b, collapse = ",")) a else b
    paste(side, collapse = "|")
  })
}
same_tree <- function(a, b) setequal(split_keys(a), split_keys(b))

# Descendant tip labels of every internal node (the rooted clades of `tree`).
clade_tipsets <- function(tree) {
  nt <- NTip(tree); labs <- TipLabels(tree)
  desc <- phangorn::Descendants(tree, (nt + 1):(nt + tree$Nnode), type = "tips")
  lapply(desc, function(ix) labs[ix])
}

run_tnt_perpass <- function(phy, seed, kpass) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  saves <- character(0)
  for (i in seq_len(kpass))
    saves <- c(saves, "sectsch=rss;", sprintf("tsave *p%d.tre;", i), "save;", "tsave/;")
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              saves, "quit;")
  writeLines(script, file.path(wd, "sharedstart.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "sharedstart.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  pass_scores <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                  grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  rd <- function(f) { t <- tryCatch(ReadTntTree(file.path(wd, f)), error = function(e) NULL)
                      if (inherits(t, "multiPhylo")) t[[1]] else t }
  list(t0 = rd("t0.tre"),
       trees = lapply(seq_len(kpass), function(i) rd(sprintf("p%d.tre", i))),
       pass_scores = pass_scores)
}

for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  n <- length(phy); smin <- round(n * 0.35); smax <- round(n * 0.65)
  for (sd in seeds) {
    cat(sprintf("\n========== %s seed %d  (n=%d, eligible clade band [%d,%d]) ==========\n",
                nm, sd, n, smin, smax))
    r <- run_tnt_perpass(phy, sd, K)
    if (is.null(r$t0)) { cat("  no T0\n"); next }
    s0 <- TreeLength(r$t0, phy)
    scores <- c(s0, r$pass_scores)
    cat(sprintf("  scores by pass: %s\n", paste(scores, collapse = " -> ")))
    # First pass that strictly drops the score
    imp <- which(diff(scores) < 0)
    if (!length(imp)) { cat("  no improving pass (sectorial found nothing)\n"); next }
    i <- imp[1]                              # 1-based pass index in r$trees
    A <- if (i == 1) r$t0 else r$trees[[i - 1]]
    B <- r$trees[[i]]
    if (is.null(A) || is.null(B)) { cat("  missing tree for pass ", i, "\n"); next }
    dropA <- scores[i] - scores[i + 1]
    cat(sprintf("  FIRST improving pass = %d : %g -> %g (drop %g)\n",
                i, scores[i], scores[i + 1], dropA))
    # RF (count of differing splits)
    kA <- split_keys(A); kB <- split_keys(B)
    rf <- length(setdiff(kA, kB)) + length(setdiff(kB, kA))
    cat(sprintf("  RF(A,B) = %d differing splits\n", rf))
    if (same_tree(A, B)) { cat("  (trees identical - score drop without topology change?!)\n"); next }
    # Single-SPR moved-set test: smallest clade C of A whose removal makes A,B match
    csets <- clade_tipsets(A)
    csets <- csets[order(lengths(csets))]
    moved <- NULL
    for (C in csets) {
      if (length(C) < 2 || length(C) > n - 2) next
      Am <- tryCatch(KeepTip(A, setdiff(TipLabels(A), C)), error = function(e) NULL)
      Bm <- tryCatch(KeepTip(B, setdiff(TipLabels(B), C)), error = function(e) NULL)
      if (!is.null(Am) && !is.null(Bm) && same_tree(Am, Bm)) { moved <- C; break }
    }
    if (is.null(moved)) {
      cat("  NOT a single clade-SPR: no clade-removal makes A==B.\n")
      cat("  => either a non-clade band, a TBR (re-rooted regraft), or multiple moves.\n")
      cat(sprintf("  lost splits (in A, |smaller side|): %s\n",
                  paste(sort(sapply(setdiff(kA, kB), function(k) length(strsplit(k,"\\|")[[1]]))), collapse=",")))
      cat(sprintf("  gained splits (in B, |smaller side|): %s\n",
                  paste(sort(sapply(setdiff(kB, kA), function(k) length(strsplit(k,"\\|")[[1]]))), collapse=",")))
    } else {
      sz <- length(moved); inband <- sz >= smin && sz <= smax
      # Attachment-only vs internal change: compare induced topology on moved set
      indA <- KeepTip(A, moved); indB <- KeepTip(B, moved)
      attach_only <- if (sz >= 4) same_tree(indA, indB) else TRUE
      cat(sprintf("  SINGLE clade-SPR. moved clade size = %d  (in-band [%d,%d]? %s)\n",
                  sz, smin, smax, inband))
      cat(sprintf("  moved-clade internal topology: %s\n",
                  if (attach_only) "UNCHANGED  => ATTACHMENT/ROOTING-only (b=1 floating-HTU fix)"
                  else "CHANGED     => internal rearrangement (RAS diversity/acceptance)"))
    }
  }
}
