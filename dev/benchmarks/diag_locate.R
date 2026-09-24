# LOCALISE THE MISSING MOVE (advisor rung 4): where do T0 (487) and TNT's best
# (482) differ on Wortley? Prune to the symmetric split-difference. If the 5 steps
# live in ONE small sub-clade, a single correct sector should recover it (=> our
# selection or reduction misses it). If spread across many clades, it's the
# ITERATION of accept-equal resolves that matters.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
nm  <- Sys.getenv("TS_DATASET", "Wortley2006")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), paste0("locate", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

phy <- fitch(inapplicable.phyData[[nm]])
WriteTntCharacters(phy, file.path(wd, "data.tnt"))
script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", "rseed 1;", "taxname=;",
            "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
            rep("sectsch=rss;", 8), "tsave *best.tre;", "save;", "tsave/;", "quit;")
writeLines(script, file.path(wd, "loctest.run"))
old <- setwd(wd)
out <- suppressWarnings(system2(TNT, args = "loctest.run;", stdout = TRUE, stderr = TRUE))
setwd(old)
t0   <- ReadTntTree(file.path(wd, "t0.tre"));   if (inherits(t0, "multiPhylo"))   t0 <- t0[[1]]
best <- ReadTntTree(file.path(wd, "best.tre")); if (inherits(best, "multiPhylo")) best <- best[[1]]
labs <- TipLabels(t0)
best <- KeepTip(best, labs)
cat(sprintf("%s: TreeLength T0=%.0f  best=%.0f  (diff %+.0f)\n",
            nm, TreeLength(t0, phy), TreeLength(best, phy),
            TreeLength(best, phy) - TreeLength(t0, phy)))

# Clades (bipartitions) of each tree as canonical tip-sets, via ape::prop.part.
clades <- function(tr) {
  pp <- ape::prop.part(tr); lab <- attr(pp, "labels")
  lapply(pp, function(ix) sort(lab[ix]))
}
small <- function(s) if (length(s) <= length(labs) / 2) s else sort(setdiff(labs, s))
key <- function(s) paste(small(s), collapse = "|")
c0 <- clades(t0); cb <- clades(best)
k0 <- vapply(c0, key, ""); kb <- vapply(cb, key, "")
gained <- cb[!(kb %in% k0)]; lost <- c0[!(k0 %in% kb)]
cat(sprintf("RF = %d  (%d clades gained, %d lost)\n",
            length(gained) + length(lost), length(gained), length(lost)))
cat("\n-- small side of each GAINED clade (the rearrangement TNT made that we lack) --\n")
involved <- character(0)
for (g in gained) { sm <- small(g); involved <- union(involved, sm)
  cat(sprintf("  [%2d tips] %s\n", length(sm), paste(sm, collapse = ", "))) }
cat(sprintf("\nUNION of tips in gained clades: %d of %d total\n", length(involved), length(labs)))
# smallest clade of T0 that CONTAINS all involved tips (the sector that must be picked)
contain_sz <- vapply(c0, function(s) { sd <- small(s)
  if (all(involved %in% sd)) length(sd) else .Machine$integer.max }, integer(1))
cat(sprintf("Smallest T0 clade containing all moved tips: %d tips (sector must cover this)\n",
            min(contain_sz)))
cat("\nT0   newick:\n"); cat(ape::write.tree(ape::ladderize(t0)),  "\n")
cat("\nbest newick:\n"); cat(ape::write.tree(ape::ladderize(best)), "\n")
