# 40,000-tree matrix-free smoke test for WideSample() over the
# MaxMin::FarFirst() distance-column oracle path.
suppressPackageStartupMessages({
  library(TreeTools)   # for as.phylo.numeric
  library(TreeSearch)
})
N <- 40000L
n <- 10L
cat("Building", N, "trees...\n")
trees <- as.phylo(0:(N - 1L), nTip = 8L)

gc(reset = TRUE)
t0 <- Sys.time()
sub <- WideSample(trees, n)            # quality = NULL -> tier 1 (matrix-free)
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
peak_mb <- sum(gc()[, 6])              # max used (Mb), Ncells + Vcells

m  <- as.matrix(TreeDist::ClusteringInfoDistance(sub))
tk <- min(m[lower.tri(m)])

cat(sprintf("RESULT: %d trees selected, Tk=%.4f, wall=%.1fs, peakRAM~%.0fMB\n",
            length(sub), tk, wall, peak_mb))
cat("(A dense 40000x40000 double matrix would be ~12800 MB; matrix-free stays far below.)\n")
