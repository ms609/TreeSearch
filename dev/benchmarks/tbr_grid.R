# tbr_grid.R -- the deliverable ensemble grid.
#
# For each dataset, build a SHARED ladder of start trees spanning a range of
# qualities, write each to Newick ONCE, then feed the IDENTICAL Newick into
# BOTH engines and run TBR-to-convergence across several seeds, in two modes:
#   Mode A  strict single-tree TBR  (TS acceptEqual=F ; TNT nomulpars hold 1)
#   Mode B  buffer / plateau        (TS acceptEqual=T ; TNT mulpars  hold 1000)
# Results are paired by start tree.  Writes a tidy CSV + prints summary tables.
#
# Usage: Rscript dev/benchmarks/tbr_grid.R [datasets...]   (default Zanol Zhu)
source("dev/benchmarks/tbr_shared_start_lib.R")

args  <- commandArgs(trailingOnly = TRUE)
DSETS <- if (length(args)) args else c("Zanol2014", "Zhu2013")
SEEDS <- 1:6
OUTDIR <- "dev/benchmarks/tbr_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

asPhylo <- function(edge, d)
  structure(list(edge = edge, Nnode = d$nTip - 1L,
                 tip.label = names(d$phy)), class = "phylo")

# Build a deterministic quality ladder of start trees for dataset `d`.
# Returns a named list of phylo trees (tips renumbered to d$phy order).
buildStartLadder <- function(d) {
  norm <- function(tr) Preorder(RenumberTips(tr, names(d$phy)))
  ladder <- list()
  # 2 random topologies (poorest)
  for (i in 1:2) {
    set.seed(1000 + i)
    ladder[[paste0("random", i)]] <- norm(RandomTree(d$phy, root = TRUE))
  }
  # 2 RAS Wagner trees (poor)
  for (i in 1:2) {
    set.seed(2000 + i)
    w <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
    ladder[[paste0("wagner", i)]] <- norm(asPhylo(w$edge, d))
  }
  # 1 partially-TBR-optimised tree (medium): a Wagner pushed ~15 accepted moves
  set.seed(3001)
  w <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
  partial <- TsTbr(d, norm(asPhylo(w$edge, d)), seed = 3001,
                   acceptEqual = FALSE, maxChanges = 15L)$tree
  ladder[["partial"]] <- norm(partial)
  # near-optimal anchor: canonical TNT T0 (both engines should hold it)
  ladder[["t0anchor"]] <- norm(ape::read.tree(file.path(T0_DIR, paste0(d$name, ".tre"))))
  ladder
}

allRows <- list()
for (dn in DSETS) {
  d <- prepareDataset(dn)
  cat(sprintf("\n=== %s (n=%d) ===\n", dn, d$nTip))
  ladder <- buildStartLadder(d)
  # Persist the shared starts as a single multi-Newick file (inspectable).
  starts <- structure(ladder, class = "multiPhylo")
  ape::write.tree(starts, file.path(OUTDIR, paste0(dn, "_starts.nwk")))
  cat("start lengths:",
      paste(sprintf("%s=%.0f", names(ladder),
                    vapply(ladder, TreeLength, double(1), d$phy)), collapse = "  "), "\n")

  for (sname in names(ladder)) {
    st <- ladder[[sname]]
    for (s in SEEDS) {
      runs <- list(
        A_strict = list(TntTbr(d, st, seed = s, mulpars = FALSE, hold = 1),
                        TsTbr(d, st, seed = s, acceptEqual = FALSE)$row),
        B_buffer = list(TntTbr(d, st, seed = s, mulpars = TRUE, hold = 1000),
                        TsTbr(d, st, seed = s, acceptEqual = TRUE, maxHits = 50L)$row))
      for (md in names(runs)) for (rr in runs[[md]]) {
        rr$dataset <- dn; rr$start <- sname; rr$mode <- md
        allRows[[length(allRows) + 1]] <- rr
      }
    }
    cat(".")
  }
  cat(" done\n")
}

res <- do.call(rbind, lapply(allRows, function(r) r[, c(
  "dataset","start","mode","engine","seed",
  "start_len","final_len","rearrangements")]))
res$improvement <- res$start_len - res$final_len
csv <- file.path(OUTDIR, paste0("tbr_grid_", paste(DSETS, collapse = "_"), ".csv"))
write.csv(res, csv, row.names = FALSE)
cat("\nWrote", nrow(res), "rows to", csv, "\n")

# ---- Summary: per dataset x start x mode x engine -> final_len distribution ----
cat("\n=== ENSEMBLE SUMMARY (final length over", length(SEEDS), "seeds) ===\n")
fmt <- function(v) sprintf("min=%.0f med=%.0f max=%.0f",
                           min(v), stats::median(v), max(v))
for (dn in unique(res$dataset)) for (md in unique(res$mode)) {
  cat(sprintf("\n-- %s  mode %s --\n", dn, md))
  sub <- res[res$dataset == dn & res$mode == md, ]
  for (sn in unique(sub$start)) {
    ss <- sub[sub$start == sn, ]
    sl <- ss$start_len[1]
    tnt <- ss$final_len[ss$engine == "TNT"]; ts <- ss$final_len[ss$engine == "TS"]
    cat(sprintf("  %-9s start=%-4.0f  TNT[%s]  TS[%s]  gap(medianTS-medianTNT)=%+.0f\n",
        sn, sl, fmt(tnt), fmt(ts), stats::median(ts) - stats::median(tnt)))
  }
}
