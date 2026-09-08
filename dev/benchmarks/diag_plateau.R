# PLATEAU EXPERIMENT (oracle, mechanism check -- NOT the gate). Does accepting
# equal-length RAS-rebuild alternatives in the sector search (sectorAcceptEqual +
# rasStarts>1) let our strict-descent sectorial escape the TBR-local optimum T0?
# Per advisor: Wortley-from-T0 dropping ANY amount below 487 = positive direction
# (Wortley already ties end-to-end). The honest gate is end-to-end Zanol < 1265,
# tested separately. SCORE is the signal (candidates_evaluated is untrustworthy).
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006")), "\\s+")[[1]]
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), paste0("plateau", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, showWarnings = FALSE, recursive = TRUE)
get_t0 <- function(phy, seed = 1) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", 8), "quit;")
  writeLines(script, file.path(wd, "plttest.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "plttest.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- ReadTntTree(file.path(wd, "t0.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, tnt = if (length(s_sect)) s_sect[length(s_sect)] else NA)
}
run <- function(d, tree, ras, ae) {
  set.seed(1)
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, rssRounds = 8L,
    rasStarts = as.integer(ras), sectorAcceptEqual = ae))
  as.double(attr(r, "score"))
}
arms <- list(
  c(1, FALSE),  # base: default polish, strict descent (control)
  c(1, TRUE),   # ae only, no rebuild -> change is inert (expect == base)
  c(3, FALSE),  # rebuild, discard equal (old behaviour, expect == base)
  c(3, TRUE),   # rebuild + keep equal alternative  <-- THE TEST
  c(10, TRUE))  # more rebuild starts -> more plateau escape routes
lab <- c("ras1/ae0", "ras1/ae1", "ras3/ae0", "ras3/ae1", "ras10/ae1")
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  tn <- get_t0(phy)
  start <- TreeLength(tn$t0, phy)
  cat(sprintf("\n==== %s | T0=%.0f  TNT_sect=%.0f ====\n", nm, start, tn$tnt))
  for (i in seq_along(arms)) {
    sc <- run(phy, tn$t0, arms[[i]][1], as.logical(arms[[i]][2]))
    cat(sprintf("  %-9s score=%.0f  (%+.0f vs T0)%s\n", lab[i], sc, sc - start,
                ifelse(sc < start, "  <-- ESCAPED", "")))
  }
}
