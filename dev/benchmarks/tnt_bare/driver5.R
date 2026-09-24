source("dev/benchmarks/tnt_bare/harness.R")
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
runbest <- function(lines) { out <- run_tnt(c(lines,"tsave *finalt.tre;","save;","tsave/;","quit;"))
  min(rss_bests(out)) }
summ <- function(v) sprintf("min=%g median=%g max=%g  {%s}", min(v), median(v), max(v), paste(v,collapse=","))
SEEDS <- 1:6

cat("==== Q1: plateau distribution over seeds (30 rounds) ====\n")
# Single fixed T0, vary sectsch seed
for (cfg in c("noglobal noequals", "equals")) {
  v <- sapply(SEEDS, function(s) runbest(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
        "hold 1000;","proc tee.tre;",sprintf("sectsch: %s;",cfg), rep("sectsch=rss;",30))))
  cat(sprintf("  SINGLE-T0  [%-18s]  %s\n", cfg, summ(v)))
}
# Fixed diverse SET (seed-1), vary sectsch seed, strict
v <- sapply(SEEDS, function(s) runbest(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
      "hold 1000;","proc set.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",30))))
cat(sprintf("  SET(10div) [%-18s]  %s\n", "noglobal noequals", v |> summ()))

cat("\n==== Q2: DIVERSITY vs EFFORT — single-T0 strict with 10x rounds (300) ====\n")
cat("     (if effort alone reached 1261, single-strict@300 ~ set-strict@30)\n")
v <- sapply(1:4, function(s) runbest(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
      "hold 1000;","proc tee.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",300))))
cat(sprintf("  SINGLE-T0 strict @300 rounds, seeds 1-4: %s\n", summ(v)))
