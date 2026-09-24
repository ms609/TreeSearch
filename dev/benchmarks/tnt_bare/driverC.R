source("dev/benchmarks/tnt_bare/harness.R")
# Mechanism test: is the set's advantage just LARGER SECTORS (shared size-increase counter
# advances faster with 10 trees), not info transfer? Force a SINGLE tree to use big sectors.
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
best <- function(lines) min(rss_bests(run_tnt(c(lines, "quit;"))))
SEEDS <- 1:6
hit <- function(v) sprintf("1261 in %d/%d, med=%g [%g-%g]", sum(v<=1261), length(v), median(v), min(v), max(v))

cat("== SINGLE tree, strict, forced sector size minsize=maxsize=K (30 rounds) ==\n")
for (K in c(37, 45, 55, 65, 70)) {
  v <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
        "hold 1000;","proc tee.tre;",
        sprintf("sectsch: noglobal noequals minsize %d maxsize %d;", K, K),
        rep("sectsch=rss;",30))))
  cat(sprintf("  size=%2d : %s\n", K, hit(v)))
}
cat("\n== references ==\n")
v <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
      "hold 1000;","proc tee.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",30))))
cat(sprintf("  SINGLE default size : %s\n", hit(v)))
v <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
      "hold 1000;","proc set.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",30))))
cat(sprintf("  SET (10 trees)      : %s\n", hit(v)))
