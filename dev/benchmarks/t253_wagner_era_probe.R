#!/usr/bin/env Rscript
# t253 DECISION PROBE: is the March->current Wagner improvement GENERAL, or is it
# project4284's alone?
#
# WHY THIS EXISTS.  `dev/benchmarks/t252_mbank_*.csv` (2026-03-27) were produced by
# an engine whose Wagner addition was several-fold worse than the current one: on
# project4284, bare `AdditionTree` with the same data and seed and NO SEARCH scores
# 1590 under the March library against 409 under the current one (worse AND 5x
# faster -- the signature of the insertion-cost bug).  `t253_conv_gap_mbank.csv` and
# `t253_gap_characterization.md` publish an n=23 Spearman analysis derived from those
# t252 cells, so the question is not academic:
#
#   * if the improvement is BROAD, every t252 start tree is systematically bad and
#     the t253 analysis rests on them -> RETRACT it.
#   * if it is essentially project4284 alone, one row is contaminated -> ANNOTATE.
#
# Nothing else decides this.  A search-level comparison would confound the addition
# tree with everything the search does afterwards, which is why this measures BARE
# `AdditionTree` and nothing else.
#
# THREE CONFOUNDS THIS DESIGN REMOVES, each of which would fake a result:
#
#  1. PREPROCESSING.  The two libraries ship different TreeTools versions, so running
#     `ReadAsPhyDat` + `PhyDatToMatrix` under each would let a preprocessing
#     difference masquerade as an engine difference.  STEP `prep` therefore builds
#     every dataset ONCE, under one library, and saves it as RDS; both engines read
#     the identical object.
#  2. THE SCORER.  Asking each engine to score its own tree compares scorers as well
#     as builders.  STEP `score` re-scores EVERY saved tree under ONE library, so the
#     numbers differ only by tree quality.  Cf. na-validation-alignment-gotcha.
#  3. GAP TREATMENT.  t252 kept "-" as a sixth level (BGS); this maps "-" -> "?"
#     (plain Fitch) for BOTH arms, matching `reach_escalation_ab.R`.  That is worth
#     +61 on one tree, so mixing the two conventions across arms would dwarf the
#     effect being measured.  Both arms use the SAME convention, which is what makes
#     this an engine comparison.
#
# Env: STEP (prep|addition|score), OUT_DIR, ENGINE (label, for `addition`),
#      NEOTRANS_DIR + CAT_CSV (for `prep`), N_SEEDS (default 3).
# The library selection is the CALLER's job, via R_LIBS -- see the .sh.

step <- Sys.getenv("STEP", "")
outDir <- Sys.getenv("OUT_DIR", "")
engine <- Sys.getenv("ENGINE", "")
nSeeds <- as.integer(Sys.getenv("N_SEEDS", "3"))
if (!nzchar(outDir)) stop("OUT_DIR unset")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
pdDir <- file.path(outDir, "pd")
treeDir <- file.path(outDir, "trees")

# The fixed 25-matrix training sample, verbatim from bench_datasets.R /
# reach_escalation_ab.R: "results are only comparable when the same sample is used".
MBANK_FIXED_SAMPLE <- c(
  "project532", "project2346", "project2451", "project4501",
  "project944", "project971_(1)", "project2762",
  "project826", "project561", "project571", "project4146_(3)",
  "project3688", "project4049", "project423",
  "project4286", "project4359", "project4397", "project2084_(1)",
  "project2771", "project2184", "project3938",
  "syab07201", "project4133", "project804", "project4284"
)

BASE_SEED <- 1301L
seedsFor <- function(i) BASE_SEED + seq_len(nSeeds) - 1L + (i - 1L) * 100L

safeKey <- function(k) gsub("[^A-Za-z0-9]", "_", k)

# ---------------------------------------------------------------- STEP: prep ----
if (identical(step, "prep")) {
  suppressPackageStartupMessages(library("TreeTools"))
  neoDir <- Sys.getenv("NEOTRANS_DIR", "")
  catCsv <- Sys.getenv("CAT_CSV", "")
  if (!nzchar(neoDir)) stop("NEOTRANS_DIR unset")
  if (!nzchar(catCsv)) stop("CAT_CSV unset")
  catalogue <- read.csv(catCsv, stringsAsFactors = FALSE)
  rownames(catalogue) <- catalogue$key
  dir.create(pdDir, showWarnings = FALSE, recursive = TRUE)

  toFitch <- function(pd) {
    m <- PhyDatToMatrix(pd, ambigNA = FALSE)
    m[m == "-"] <- "?"
    MatrixToPhyDat(m)
  }

  rows <- list()
  for (key in MBANK_FIXED_SAMPLE) {
    if (!key %in% catalogue$key) {
      message("SKIP (not in catalogue): ", key)
      next
    }
    row <- catalogue[key, ]
    # SEQUESTER the validation split: a one-way door.
    if (!identical(row$split, "training")) {
      message("SKIP (split=", row$split, ", SEQUESTERED): ", key)
      next
    }
    f <- file.path(neoDir, row$filename)
    if (!file.exists(f)) {
      message("SKIP (file missing): ", key)
      next
    }
    pd <- toFitch(suppressWarnings(ReadAsPhyDat(f)))
    saveRDS(pd, file.path(pdDir, paste0(safeKey(key), ".rds")))
    rows[[length(rows) + 1L]] <- data.frame(
      key = key, nTip = length(pd), nChar = attr(pd, "nr"),
      stringsAsFactors = FALSE
    )
    cat(sprintf("prep %-18s nTip=%5d nChar=%4d\n", key, length(pd), attr(pd, "nr")))
  }
  write.csv(do.call(rbind, rows), file.path(outDir, "manifest.csv"),
            row.names = FALSE)
  cat("prep done:", length(rows), "datasets\n")
}

# ------------------------------------------------------------ STEP: addition ----
# Bare `AdditionTree`, no search.  Trees are saved as RDS rather than Newick: a
# round-trip through Newick could renumber or reorder, and STEP `score` must see the
# tree the engine actually produced.
if (identical(step, "addition")) {
  if (!nzchar(engine)) stop("ENGINE unset")
  suppressPackageStartupMessages(library("TreeSearch"))
  dir.create(treeDir, showWarnings = FALSE, recursive = TRUE)
  tsVer <- as.character(utils::packageVersion("TreeSearch"))
  tsLib <- dirname(system.file(package = "TreeSearch"))
  cat("engine=", engine, " TreeSearch ", tsVer, " from ", tsLib, "\n", sep = "")

  manifest <- read.csv(file.path(outDir, "manifest.csv"), stringsAsFactors = FALSE)
  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    key <- manifest$key[i]
    pd <- readRDS(file.path(pdDir, paste0(safeKey(key), ".rds")))
    for (sd in seedsFor(i)) {
      set.seed(sd)
      t0 <- Sys.time()
      tr <- tryCatch(AdditionTree(pd), error = function(e) e)
      wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
      if (inherits(tr, "error")) {
        cat(sprintf("FAIL %-18s seed=%d: %s\n", key, sd, conditionMessage(tr)))
        rows[[length(rows) + 1L]] <- data.frame(
          key = key, engine = engine, tsVersion = tsVer, tsLib = tsLib, seed = sd,
          wallS = wall, ok = FALSE, stringsAsFactors = FALSE
        )
        next
      }
      # Write the deliverable BEFORE anything that could fail on it: a verification
      # error must not take the computed tree with it (the 2026-07-31 lesson).
      saveRDS(tr, file.path(treeDir, sprintf("%s__%s__s%d.rds",
                                             safeKey(key), engine, sd)))
      rows[[length(rows) + 1L]] <- data.frame(
        key = key, engine = engine, tsVersion = tsVer, tsLib = tsLib, seed = sd,
        wallS = wall, ok = TRUE, stringsAsFactors = FALSE
      )
      cat(sprintf("add  %-18s %-6s seed=%d  %.2fs\n", key, engine, sd, wall))
    }
  }
  write.csv(do.call(rbind, rows),
            file.path(outDir, paste0("addition_", engine, ".csv")),
            row.names = FALSE)
}

# --------------------------------------------------------------- STEP: score ----
# ONE scorer for every tree, whichever engine built it.
if (identical(step, "score")) {
  suppressPackageStartupMessages(library("TreeSearch"))
  cat("scorer: TreeSearch ", as.character(utils::packageVersion("TreeSearch")),
      "\n", sep = "")
  manifest <- read.csv(file.path(outDir, "manifest.csv"), stringsAsFactors = FALSE)

  # ---- PROVENANCE ASSERTION: both engines report Version 2.0.0, so the version
  # string CANNOT distinguish them.  If R_LIBS failed to take, both arms would run
  # the SAME library and every ratio would be 1.000 -- which reads exactly like the
  # "annotate, not retract" answer.  A null result must not be obtainable from a
  # broken arm, so refuse to report unless the two arms provably used different
  # library paths.
  addFiles <- list.files(outDir, pattern = "^addition_.*\\.csv$", full.names = TRUE)
  if (length(addFiles) >= 2L) {
    prov <- do.call(rbind, lapply(addFiles, read.csv, stringsAsFactors = FALSE))
    byEngine <- unique(prov[, c("engine", "tsVersion", "tsLib")])
    cat("\n---- provenance ----\n")
    print(byEngine, row.names = FALSE)
    if (anyDuplicated(byEngine$tsLib)) {
      stop("PROVENANCE FAILURE: two engine labels resolved to the SAME library ",
           "path, so the arms are not independent. R_LIBS did not take. Any ratio ",
           "of 1.000 here would be an artefact, not a finding.")
    }
    cat("provenance OK: arms used distinct libraries\n")
  }

  files <- list.files(treeDir, pattern = "\\.rds$", full.names = TRUE)
  rows <- list()
  for (f in files) {
    parts <- strsplit(sub("\\.rds$", "", basename(f)), "__", fixed = TRUE)[[1]]
    if (length(parts) != 3L) {
      cat("SKIP unparseable filename:", basename(f), "\n")
      next
    }
    keySafe <- parts[[1]]
    idx <- match(keySafe, safeKey(manifest$key))
    if (is.na(idx)) {
      cat("SKIP no manifest row:", basename(f), "\n")
      next
    }
    key <- manifest$key[idx]
    pd <- readRDS(file.path(pdDir, paste0(keySafe, ".rds")))
    tr <- readRDS(f)
    sc <- tryCatch(TreeLength(tr, pd, concavity = Inf),
                   error = function(e) {
                     cat("SCORE FAIL", basename(f), ":", conditionMessage(e), "\n")
                     NA_real_
                   })
    rows[[length(rows) + 1L]] <- data.frame(
      key = key, engine = parts[[2]], seed = as.integer(sub("^s", "", parts[[3]])),
      nTip = manifest$nTip[idx], nChar = manifest$nChar[idx],
      score = as.double(sc), stringsAsFactors = FALSE
    )
  }
  scores <- do.call(rbind, rows)
  write.csv(scores, file.path(outDir, "scores.csv"), row.names = FALSE)

  # ------- the decision table: per MATRIX, not per (matrix, seed) -------
  # Pairing on cells would be pseudo-replicated -- seeds within a matrix are not
  # independent evidence about the ENGINE.
  # ORIENTATION IS LOAD-BEARING.  The ratio must be OLD / NEW, so that "> 1" means
  # "the March engine built a worse (longer) tree" -- which is what every sentence
  # below asserts.  Alphabetical order would put `cur` in the numerator and silently
  # invert every count against its own caption.
  engines <- unique(scores$engine)
  engines <- c(intersect(c("t252"), engines), sort(setdiff(engines, "t252")))
  if (length(engines) == 2L) {
    medOf <- function(k, e) {
      v <- scores$score[scores$key == k & scores$engine == e]
      if (!length(v) || all(is.na(v))) NA_real_ else median(v, na.rm = TRUE)
    }
    keys <- unique(scores$key)
    tab <- data.frame(
      key = keys,
      nTip = manifest$nTip[match(keys, manifest$key)],
      a = vapply(keys, medOf, double(1), e = engines[[1]]),
      b = vapply(keys, medOf, double(1), e = engines[[2]]),
      stringsAsFactors = FALSE
    )
    names(tab)[3:4] <- engines
    tab$ratio <- tab[[3]] / tab[[4]]
    tab <- tab[order(-tab$ratio), ]
    write.csv(tab, file.path(outDir, "decision_table.csv"), row.names = FALSE)
    cat("\n==== per-MATRIX medians (", engines[[1]], " vs ", engines[[2]],
        ") ====\n", sep = "")
    print(tab, row.names = FALSE)
    ok <- !is.na(tab$ratio)
    cat("\nmatrices where", engines[[1]], "is WORSE (ratio > 1):",
        sum(tab$ratio[ok] > 1), "of", sum(ok), "\n")
    cat("matrices within 1%:", sum(abs(tab$ratio[ok] - 1) < 0.01), "\n")
    cat("matrices >10% worse:", sum(tab$ratio[ok] > 1.10), "\n")
    cat("median ratio:", median(tab$ratio[ok]), "\n")
    cat("\nREAD THIS AS: broad ratios > 1 => the t253 n=23 analysis rests on",
        "systematically bad start trees (RETRACT).  A ratio > 1 on project4284",
        "alone => ANNOTATE that row.\n")
  } else {
    cat("\nOnly", length(engines), "engine(s) present; run both arms.\n")
  }
}
