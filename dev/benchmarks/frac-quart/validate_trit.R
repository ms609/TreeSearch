# Validation of the unit = "nrqs" path in QuartetConcordance().
#
# Ordering follows the advisor's priority:
#   1. Cell extraction cross-checked against the verified C++ kernel (conc/dec).
#   2. A_ij <= min(Wc, Wk) asserted on every real (split, char, pair).
#   3. Spec ladder (t = 8) reproduced from the per-pair formula.
#   4. An INDEPENDENT re-implementation of the measure matches the package on
#      binary (congreveLamsdell), multistate, and missing-data inputs.
#   5. `return` parsing quirk mirrored between quartet and trit paths.
#   6. congreveLamsdell quartet -> trit edge-score movement + regime census.
#
# Run: Rscript dev/benchmarks/frac-quart/validate_trit.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
library("TreeTools", quietly = TRUE)

ok <- function(label, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) stop("FAILED: ", label)
}

# ---- Rebuild the (logiSplits, charInt) the function feeds to its kernels ----
build_inputs <- function(tree, dataset) {
  tipLabels <- intersect(TipLabels(tree), names(dataset))
  dataset <- dataset[tipLabels, drop = FALSE]
  splits <- as.Splits(tree, dataset)
  logiSplits <- vapply(seq_along(splits), function(i) as.logical(splits[[i]]),
                       logical(NTip(dataset)))
  contrast <- attr(dataset, "contrast")
  charLevels <- attr(dataset, "allLevels")
  isInapp <- charLevels == "-"
  isAmbig <- rowSums(contrast[, colnames(contrast) != "-", drop = FALSE]) > 1
  isGrouping <- !isAmbig & !isInapp
  groupingCols <- apply(contrast[isGrouping, , drop = FALSE] > 0, 1, which)
  levelToInt <- rep(NA_integer_, length(charLevels))
  levelToInt[isGrouping] <- as.integer(groupingCols)
  characters <- PhyDatToMatrix(dataset)
  charInt <- array(levelToInt[match(characters, charLevels)],
                   dim = dim(characters), dimnames = dimnames(characters))
  list(logiSplits = logiSplits, charInt = charInt, splits = splits)
}

ch2 <- function(x) x * (x - 1) / 2
posm <- function(z) { z[z < 0] <- 0; z }

# Per (split, char) cells summed over state-pairs, returned as component lists
# so the same crosstab feeds both the kernel cross-check and the trit measure.
pair_cells <- function(logiSplits, charInt) {
  nSplit <- ncol(logiSplits); nChar <- ncol(charInt)
  out <- vector("list", nChar)
  for (ci in seq_len(nChar)) {
    col <- charInt[, ci]; scored <- !is.na(col)
    states <- sort(unique(col[scored]))
    prs <- list()
    if (length(states) >= 2) {
      for (a in seq_len(length(states) - 1L)) for (b in seq(a + 1L, length(states))) {
        inI <- scored & col == states[a]; inJ <- scored & col == states[b]
        prs[[length(prs) + 1L]] <- list(
          aI = colSums(logiSplits & inI), bI = colSums(!logiSplits & inI),
          aJ = colSums(logiSplits & inJ), bJ = colSums(!logiSplits & inJ))
      }
    }
    out[[ci]] <- prs
  }
  out
}

# ---------- 1. Cell extraction vs the verified kernel ----------
{
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  tree <- TreeSearch::referenceTree
  inp <- build_inputs(tree, dataset)
  cells <- pair_cells(inp$logiSplits, inp$charInt)
  nSplit <- ncol(inp$logiSplits); nChar <- ncol(inp$charInt)
  reconC <- reconD <- matrix(0, nSplit, nChar)
  for (ci in seq_len(nChar)) for (p in cells[[ci]]) {
    reconC[, ci] <- reconC[, ci] + ch2(p$aI) * ch2(p$bJ) + ch2(p$bI) * ch2(p$aJ)
    reconD[, ci] <- reconD[, ci] + p$aI * p$bI * p$aJ * p$bJ
  }
  reconD <- reconD + reconC
  kern <- TreeSearch:::quartet_concordance(inp$logiSplits, inp$charInt)
  ok("cell -> concordant matches kernel", all(abs(reconC - kern$concordant) < 1e-9))
  ok("cell -> decisive   matches kernel", all(abs(reconD - kern$decisive) < 1e-9))

  # ---------- 2. A_ij <= min(Wc, Wk) on every pair ----------
  worst <- -Inf
  for (ci in seq_len(nChar)) for (p in cells[[ci]]) {
    aij <- posm(p$aI - 1) * posm(p$bJ - 1) + posm(p$bI - 1) * posm(p$aJ - 1)
    nI <- p$aI + p$bI; nJ <- p$aJ + p$bJ; mA <- p$aI + p$aJ; tP <- nI + nJ
    wc <- posm(nI - 1) * posm(nJ - 1); wk <- posm(mA - 1) * posm(tP - mA - 1)
    worst <- max(worst, max(aij - pmin(wc, wk)))
  }
  ok("A_ij <= min(Wc, Wk) everywhere", worst <= 1e-9)
}

# ---------- 3. Spec ladder (t = 8) from the per-pair formula ----------
pairQ <- function(p, q, r, s) {           # edge quality A / Wk for one 2x2 table
  A <- posm(p - 1) * posm(s - 1) + posm(q - 1) * posm(r - 1)
  m <- p + r; t <- p + q + r + s
  Wk <- posm(m - 1) * posm(t - m - 1)
  A / Wk
}
ladder <- rbind(
  identical    = c(4, 0, 0, 4),
  nested       = c(3, 0, 2, 3),
  mildCrossing = c(3, 1, 1, 3),
  maxCrossing  = c(2, 2, 2, 2))
Qlad <- apply(ladder, 1, function(x) pairQ(x[1], x[2], x[3], x[4]))
expected <- c(identical = 1, nested = 0.5, mildCrossing = 4 / 9, maxCrossing = 2 / 9)
ok("ladder reproduces spec table", all(abs(Qlad - expected) < 1e-9))
ok("only identical scores 1", sum(abs(Qlad - 1) < 1e-9) == 1L)
ok("crossing < nested", Qlad["mildCrossing"] < Qlad["nested"] &&
     Qlad["maxCrossing"] < Qlad["mildCrossing"])

# ---------- 4. Independent re-implementation vs the package ----------
# Deliberately different code shape: per-pair matrices, explicit pooling.
ref_trit <- function(tree, dataset, weight = TRUE, return = "edge") {
  inp <- build_inputs(tree, dataset)
  cells <- pair_cells(inp$logiSplits, inp$charInt)
  nSplit <- ncol(inp$logiSplits); nChar <- ncol(inp$charInt)
  numE <- numC <- den <- matrix(0, nSplit, nChar); wcTot <- numeric(nChar)
  for (ci in seq_len(nChar)) for (p in cells[[ci]]) {
    aij <- posm(p$aI - 1) * posm(p$bJ - 1) + posm(p$bI - 1) * posm(p$aJ - 1)
    nI <- p$aI + p$bI; nJ <- p$aJ + p$bJ; mA <- p$aI + p$aJ; tP <- nI + nJ
    wc <- posm(nI - 1) * posm(nJ - 1); wk <- posm(mA - 1) * posm(tP - mA - 1)
    m <- pmin(wc, wk)
    den[, ci] <- den[, ci] + m
    numE[, ci] <- numE[, ci] + ifelse(wk > 0, m * aij / wk, 0)
    numC[, ci] <- numC[, ci] + ifelse(wc > 0, m * aij / wc, 0)
    wcTot[ci] <- wcTot[ci] + wc[1]
  }
  info <- wcTot > 0
  ret <- pmatch(tolower(trimws(return)), c("character", "site", "default"),
                nomatch = 3L)
  if (ret == 3L) {                                  # edge
    if (isTRUE(weight)) {
      d <- rowSums(den); v <- ifelse(d == 0, NA_real_, rowSums(numE) / d)
    } else {
      sE <- ifelse(den > 0, numE / den, NA_real_)
      v <- if (any(info)) rowMeans(sE[, info, drop = FALSE], na.rm = TRUE)
           else rep(NA_real_, nSplit)
      v[is.nan(v)] <- NA_real_
    }
    setNames(v, names(inp$splits))
  } else {                                          # char
    if (isTRUE(weight)) {
      d <- colSums(den); ifelse(d == 0, NA_real_, colSums(numC) / d)
    } else {
      sC <- ifelse(den > 0, numC / den, NA_real_)
      vapply(seq_len(nChar), function(ci) if (info[ci]) {
        mm <- mean(sC[, ci], na.rm = TRUE); if (is.nan(mm)) NA_real_ else mm
      } else NA_real_, double(1))
    }
  }
}

same <- function(a, b) all((is.na(a) & is.na(b)) | (abs(a - b) < 1e-9), na.rm = FALSE)

data("congreveLamsdellMatrices", package = "TreeSearch")
binDat <- congreveLamsdellMatrices[[1]]
refTree <- TreeSearch::referenceTree
for (w in c(TRUE, FALSE)) for (r in c("edge", "char")) {
  pkg <- QuartetConcordance(refTree, binDat, weight = w, return = r, unit = "nrqs", chanceCorrect = FALSE)
  rf  <- ref_trit(refTree, binDat, weight = w, return = r)
  ok(sprintf("binary  pkg==ref  weight=%s return=%s", w, r), same(pkg, rf))
}

# Multistate toy (3- and 4-state characters, complete)
msMat <- matrix(c(0, 0, 1, 1, 2, 2, 2, 0,     # 3 states
                  0, 0, 0, 1, 1, 2, 3, 3,     # 4 states
                  0, 0, 1, 1, 2, 2, 0, 1), 8, # 3 states
                dimnames = list(paste0("t", 1:8), NULL))
msDat <- MatrixToPhyDat(msMat)
msTree <- BalancedTree(8)
for (w in c(TRUE, FALSE)) for (r in c("edge", "char")) {
  pkg <- QuartetConcordance(msTree, msDat, weight = w, return = r, unit = "nrqs", chanceCorrect = FALSE)
  rf  <- ref_trit(msTree, msDat, weight = w, return = r)
  ok(sprintf("multi   pkg==ref  weight=%s return=%s", w, r), same(pkg, rf))
}

# Missing / ambiguous / inapplicable toy
naMat <- matrix(c(0, 0, 1, 1, "?", "?", 0, 1,
                  0, 1, "-", 1, 0, "(01)", 1, 0,
                  0, 0, 0, 1, 1, 2, "?", 2), 8,
                dimnames = list(paste0("t", 1:8), NULL))
naDat <- MatrixToPhyDat(naMat)
for (w in c(TRUE, FALSE)) for (r in c("edge", "char")) {
  pkg <- QuartetConcordance(msTree, naDat, weight = w, return = r, unit = "nrqs", chanceCorrect = FALSE)
  rf  <- ref_trit(msTree, naDat, weight = w, return = r)
  ok(sprintf("missing pkg==ref  weight=%s return=%s", w, r), same(pkg, rf))
}

# ---------- 5. `return` parsing mirrors the quartet path ----------
e1 <- QuartetConcordance(refTree, binDat, return = "edge",      unit = "nrqs", chanceCorrect = FALSE)
e2 <- QuartetConcordance(refTree, binDat, return = "default",   unit = "nrqs", chanceCorrect = FALSE)
c1 <- QuartetConcordance(refTree, binDat, return = "char",      unit = "nrqs", chanceCorrect = FALSE)
c2 <- QuartetConcordance(refTree, binDat, return = "character", unit = "nrqs", chanceCorrect = FALSE)
c3 <- QuartetConcordance(refTree, binDat, return = "site",      unit = "nrqs", chanceCorrect = FALSE)
ok("return edge == default", same(e1, e2))
ok("return char == character == site", same(c1, c2) && same(c1, c3))
ok("edge length == n splits, char length == n chars",
   length(e1) == length(inp <- as.Splits(refTree)) &&
     length(c1) == sum(attr(binDat, "weight")))

# ---------- 6. congreveLamsdell quartet -> trit movement + census ----------
qEdge <- QuartetConcordance(refTree, binDat, unit = "quartet", chanceCorrect = FALSE)
tEdge <- QuartetConcordance(refTree, binDat, unit = "nrqs", chanceCorrect = FALSE)
cat("\nEdge concordance, quartet vs trit (congreveLamsdell[[1]]):\n")
print(round(rbind(quartet = qEdge, trit = tEdge, delta = tEdge - qEdge), 3))
cat(sprintf("\nmean quartet = %.3f  mean trit = %.3f  (trit is stricter)\n",
            mean(qEdge, na.rm = TRUE), mean(tEdge, na.rm = TRUE)))

# Regime census per edge: identical / nested / crossing state-pairs
inpc <- build_inputs(refTree, binDat)
cells <- pair_cells(inpc$logiSplits, inpc$charInt)
regime <- matrix(0L, ncol(inpc$logiSplits), 3,
                 dimnames = list(names(inpc$splits),
                                 c("identical", "nested", "crossing")))
for (ci in seq_along(cells)) for (p in cells[[ci]]) {
  empt <- (p$aI == 0) + (p$bI == 0) + (p$aJ == 0) + (p$bJ == 0)
  regime[, "identical"] <- regime[, "identical"] + (empt == 2)
  regime[, "nested"]    <- regime[, "nested"]    + (empt == 1)
  regime[, "crossing"]  <- regime[, "crossing"]  + (empt == 0)
}
cat("\nState-pair regime census per edge:\n")
print(regime)

cat("\nAll checks passed.\n")
