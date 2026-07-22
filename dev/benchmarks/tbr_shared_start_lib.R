# tbr_shared_start_lib.R
#
# Shared helpers for the isolated-TBR head-to-head between TreeSearch and
# TNT 1.6 from IDENTICAL starting trees.  Loaded by the pilot and the full
# grid driver.  See dev/plans/2026-06-18-tbr-shared-start.md for the design.
#
# Both engines optimise the SAME Fitch objective because the matrices have
# inapplicable tokens replaced by '?'.  Lengths are therefore directly
# comparable (TreeLength vs TNT `length`).

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})

TNT_EXE <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
T0_DIR  <- Sys.getenv("T0_DIR", "C:/Users/pjjg18/GitHub/TreeSearch/dev/benchmarks/t0")

# ---------------------------------------------------------------------------
# Dataset preparation: phyDat -> the array bundle ts_tbr_diagnostics wants.
# ---------------------------------------------------------------------------
prepareDataset <- function(name) {
  phy <- readRDS(file.path(T0_DIR, paste0(name, ".phy.rds")))
  at  <- attributes(phy)
  list(
    name     = name,
    phy      = phy,
    contrast = at$contrast,
    tip_data = matrix(unlist(phy, use.names = FALSE),
                      nrow = length(phy), byrow = TRUE),
    weight   = at$weight,
    levels   = at$levels,
    nTip     = length(phy)
  )
}

# ---------------------------------------------------------------------------
# TNT helpers
# ---------------------------------------------------------------------------

# ape Newick -> TNT parenthetical (space-separated, no branch lengths,
# no trailing semicolon).
ToTntTree <- function(tr) {
  nw <- ape::write.tree(tr)
  nw <- gsub(";", "", nw, fixed = TRUE)   # drop trailing ';'
  nw <- gsub(",", " ", nw, fixed = TRUE)  # commas -> spaces
  nw
}

# Run a TNT script (character vector of lines) in a fresh temp dir that
# already contains a data.tnt for `phy`.  Returns sanitised stdout lines.
# `files` is a named list of extra files to write into the working dir
# (name = filename, value = character vector of lines).
RunTnt <- function(phy, scriptLines, tag = "tnt", files = list()) {
  wd <- file.path(tempdir(), paste0(tag, Sys.getpid()))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  for (nm in names(files)) writeLines(files[[nm]], file.path(wd, nm))
  writeLines(scriptLines, file.path(wd, "swapper.run"))
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(TNT_EXE, args = "swapper.run;",
                                  stdout = TRUE, stderr = TRUE))
  iconv(out, from = "", to = "UTF-8", sub = "")
}

# Pull a single number from the first line matching `pat` (with one capture
# group), stripping thousands separators.
GrepNum <- function(out, pat) {
  hit <- grep(pat, out, value = TRUE)
  if (!length(hit)) return(NA_real_)
  suppressWarnings(as.numeric(gsub(",", "", sub(pat, "\\1", hit[1]))))
}

# Run TNT TBR (bbreak) from `startTree`, save result, read it back, score in
# R with TreeLength.  Returns a one-row data.frame.
#   mulpars/hold : equal-tree buffer controls (Mode A: FALSE/1; Mode B: TRUE/1000)
#   randclip     : randomise clip order using rseed (the stochasticity knob)
TntTbr <- function(d, startTree, seed, mulpars, hold, randclip = TRUE) {
  swap <- paste0("bbreak = tbr ",
                 if (randclip) "randclip " else "norandclip ",
                 if (mulpars) "mulpars" else "nomulpars", ";")
  script <- c("mxram 1024;", "taxname=;", "proc data.tnt;",
              paste0("rseed ", seed, ";"),
              paste0("hold ", hold, ";"),
              paste0("tread ", ToTntTree(startTree), ";"),
              swap,
              "tsave *out.tre;", "save;", "tsave/;",
              "quit;")
  wd <- file.path(tempdir(), paste0("tnttbr", Sys.getpid()))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(d$phy, file.path(wd, "data.tnt"))
  writeLines(script, file.path(wd, "swapper.run"))
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  # Wall = TNT process time only (data.tnt already written; ReadTntTree excluded).
  # Includes TNT startup/proc/tread, which INFLATES TNT wall ⇒ tnt rate is a
  # conservative LOWER bound (works in TNT's disfavour for the throughput read).
  .t0 <- Sys.time()
  out <- suppressWarnings(system2(TNT_EXE, args = "swapper.run;",
                                  stdout = TRUE, stderr = TRUE))
  .wall <- as.double(difftime(Sys.time(), .t0, units = "secs"))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")

  startScore  <- GrepNum(out, ".*Start swapping from .* \\(score ([0-9]+)\\).*")
  bestStdout  <- GrepNum(out, ".*Best score \\(TBR\\):\\s*([0-9]+).*")
  rearr       <- GrepNum(out, ".*Total rearrangements examined:\\s*([0-9,]+).*")
  # Authoritative final score: read saved tree(s), score in R (identical engine)
  trees <- tryCatch(ReadTntTree(file.path(wd, "out.tre")), error = function(e) NULL)
  finalR <- if (is.null(trees)) NA_real_ else {
    if (inherits(trees, "multiPhylo"))
      min(vapply(trees, TreeLength, double(1), d$phy)) else TreeLength(trees, d$phy)
  }
  nTrees <- if (is.null(trees)) NA_integer_ else
            if (inherits(trees, "multiPhylo")) length(trees) else 1L
  bestTree <- if (is.null(trees)) NULL else if (inherits(trees, "multiPhylo")) {
    trees[[which.min(vapply(trees, TreeLength, double(1), d$phy))]]
  } else trees
  row <- data.frame(engine = "TNT", seed = seed, mulpars = mulpars, hold = hold,
             start_len = TreeLength(startTree, d$phy),
             start_len_tnt = startScore, final_len = finalR,
             final_len_tnt = bestStdout, n_trees = nTrees,
             rearrangements = rearr, wall = .wall, stringsAsFactors = FALSE)
  attr(row, "tree") <- bestTree
  row
}

# ---------------------------------------------------------------------------
# TreeSearch helpers
# ---------------------------------------------------------------------------

# phylo -> kernel edge matrix (standard ape numbering, tips matching d$phy).
PhyloToKernelEdge <- function(tree, d) {
  tree <- RenumberTips(tree, names(d$phy))
  tree <- Preorder(tree)
  tree[["edge"]]
}

# Run TreeSearch TBR to convergence from `startTree`.  acceptEqual=FALSE is
# strict descent (Mode A); acceptEqual=TRUE plateau-walks the single tree
# (Mode B analogue).  Returns a one-row data.frame plus the pass trajectory.
TsTbr <- function(d, startTree, seed, acceptEqual, maxHits = 1L, maxChanges = 0L) {
  edge <- PhyloToKernelEdge(startTree, d)
  set.seed(seed)
  .t0 <- Sys.time()
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = maxHits, acceptEqual = acceptEqual, maxChanges = maxChanges)
  .wall <- as.double(difftime(Sys.time(), .t0, units = "secs"))
  resTree <- structure(list(edge = res$edge, Nnode = d$nTip - 1L,
                            tip.label = names(d$phy)), class = "phylo")
  finalR <- TreeLength(resTree, d$phy)
  row <- data.frame(engine = "TS", seed = seed, mulpars = NA, hold = NA,
                     start_len = TreeLength(startTree, d$phy),
                     start_len_tnt = NA_real_, final_len = finalR,
                     final_len_tnt = res$score, n_trees = 1L,
                     rearrangements = res$n_evaluated, wall = .wall,
                     stringsAsFactors = FALSE)
  attr(row, "tree") <- resTree
  list(row = row, tree = resTree,
       passes = res$passes, n_accepted = res$n_accepted, converged = res$converged)
}
