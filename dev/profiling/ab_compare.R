# Compare base vs step1 A/B runs: behaviour-neutrality (score + cand identical)
# and end-to-end wall delta.  Reads dev/profiling/ab_latest.csv (two tags).
d <- read.csv("dev/profiling/ab_latest.csv", stringsAsFactors = FALSE)
w <- reshape(d, idvar = c("dataset", "obj", "seed"), timevar = "tag",
             direction = "wide")
w$score_ok <- w$score.base == w$score.step1
w$cand_ok  <- w$cand.base  == w$cand.step1
cat("=== Behaviour-neutrality gate ===\n")
cat("score identical across all runs:", all(w$score_ok), "\n")
cat("candidates_evaluated identical:  ", all(w$cand_ok), "\n")
if (!all(w$score_ok & w$cand_ok)) {
  cat("\n!!! MISMATCH (change is NOT behaviour-neutral) !!!\n")
  print(w[!(w$score_ok & w$cand_ok),
          c("dataset","obj","seed","score.base","score.step1","cand.base","cand.step1")],
        row.names = FALSE)
}
cat("\n=== End-to-end wall A/B (sum over seeds; ratio = step1/base) ===\n")
agg <- aggregate(cbind(wall.base, wall.step1) ~ dataset + obj, w, sum)
agg$ratio <- round(agg$wall.step1 / agg$wall.base, 3)
print(agg, row.names = FALSE)
cat(sprintf("\nTOTAL wall  base %.1f s  step1 %.1f s  ratio %.3f\n",
            sum(w$wall.base), sum(w$wall.step1),
            sum(w$wall.step1) / sum(w$wall.base)))
