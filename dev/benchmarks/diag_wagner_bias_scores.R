# Does the insertion-cost bug also hit the BIASED Wagner variants (the
# production "default" strategy uses wagnerBias=1 Goloboff, temp 0.3)?
# Bias only changes the addition ORDER, not the cost formula, so we expect
# all three to be ~+350 over optimum if the formula is the culprit.
suppressMessages({
  library(TreeSearch, lib.loc = "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem")
  library(TreeTools)
})
nm  <- Sys.getenv("DS", "Zanol2014")
phy <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
n   <- length(phy)
at <- attributes(phy); contrast <- at$contrast
tipData <- matrix(unlist(phy, use.names = FALSE), nrow = n, byrow = TRUE)
weight  <- TreeSearch:::.ScaleWeight(at$weight); levels <- at$levels

set.seed(1)
biasName <- c("RANDOM", "GOLOBOFF(default)", "ENTROPY")
cat(sprintf("== %s | Wagner score by bias (n_reps=12, no TBR) ==\n", nm))
cat("(exact-insertion RAS ~1300; optimum ~1261)\n")
for (b in 0:2) {
  temp <- if (b == 0) 1.0 else 0.3
  res <- TreeSearch:::ts_wagner_bias_bench(
    contrast, tipData, weight, levels, integer(0), -1.0,
    b, temp, 12L, FALSE)
  s <- res$wagner_score
  cat(sprintf("  %-18s mean=%.0f sd=%.0f min=%.0f max=%.0f\n",
              biasName[b + 1], mean(s), sd(s), min(s), max(s)))
}
