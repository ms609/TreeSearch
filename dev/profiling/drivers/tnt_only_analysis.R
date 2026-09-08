# Where do TNT's extra MPTs come from?  Step 1 (foundational): regenerate TS's
# 256 fresh, compare to RAW zhu_tnt.tre (NOT the overwritten pooled rds), and
# score the TNT-only trees with TS's own scorer (by label).
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(TreeTools))
base <- "C:/Users/pjjg18/GitHub/wide-sample/dev/zhu2013"
dat0 <- ReadAsPhyDat(file.path(base, "zhu2013_orig.nex"))
taxa <- names(dat0)
datF <- TreeSearch:::.GapsAsMissing(dat0)          # exactly what the search used
og <- c("Galeaspida", "Osteostraci")
prep <- function(trs) {
  if (inherits(trs, "phylo")) trs <- c(trs)
  trs <- lapply(trs, function(t) SortTree(RootTree(t, og)))
  class(trs) <- "multiPhylo"; trs
}

# --- TNT trees (raw) ---
tnt <- prep(ReadTntTree(file.path(base, "zhu_tnt.tre"), tipLabels = taxa))

# --- TS 256 (fresh, the exhaustion-confirmed settings) ---
set.seed(2013L)
ts <- MaximizeParsimony(dat0, concavity = Inf, inapplicable = "missing",
  strategy = "intensive", maxReplicates = 15L, maxSeconds = 150,
  control = SearchControl(poolMaxSize = 20000L, tbrMaxHits = 50L,
                          enumTimeFraction = 0.5),
  verbosity = 0L)
tsP <- prep(ts)

kt <- vapply(tnt, ape::write.tree, character(1))
ks <- vapply(tsP, ape::write.tree, character(1))
tnt_only <- tnt[!(kt %in% ks)]; class(tnt_only) <- "multiPhylo"
ts_only  <- tsP[!(ks %in% kt)]; class(ts_only)  <- "multiPhylo"
cat(sprintf("TNT %d | TS %d | shared %d | TNT-only %d | TS-only %d | union %d\n",
            length(tnt), length(tsP), sum(kt %in% ks),
            length(tnt_only), length(ts_only),
            length(unique(c(kt, ks)))))

# --- score TNT-only by TS's scorer (by label) ---
sc_tnt_only <- vapply(tnt_only, function(t) TreeLength(t, datF, concavity = Inf), numeric(1))
sc_ts <- vapply(tsP[1:5], function(t) TreeLength(t, datF, concavity = Inf), numeric(1))
cat("TS pool sample scores:", paste(sc_ts, collapse = ","), "\n")
cat("TNT-only scores: range", min(sc_tnt_only), "-", max(sc_tnt_only),
    "| all == 598:", all(sc_tnt_only == 598), "\n")
print(table(sc_tnt_only))
saveRDS(list(tnt = tnt, ts = tsP, tnt_only = tnt_only, ts_only = ts_only,
             sc_tnt_only = sc_tnt_only),
        "dev/profiling/zhu_tnt_ts_compare.rds")
cat("saved dev/profiling/zhu_tnt_ts_compare.rds\n")
