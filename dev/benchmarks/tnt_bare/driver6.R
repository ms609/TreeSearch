source("dev/benchmarks/tnt_bare/harness.R")
# Lucky-tree control: does ANY single tree of the 10-tree 1271 set reach target SOLO under
# strict sectsch? If not, the set's escape needs cross-set variety, not one lucky member.
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
trees <- grep("^\\(", L, value = TRUE)              # the 10 newick lines (some end '*', last ';')
trees <- sub("[*;]$", "", trees)
runbest <- function(lines) { out <- run_tnt(c(lines,"quit;")); min(rss_bests(out)) }

cat(sprintf("==== %s: each of %d set trees, SOLO single-tree strict sectsch (seeds 1-3, 30 rounds) ====\n",
            nm, length(trees)))
solo <- numeric(length(trees))
for (i in seq_along(trees)) {
  writeLines(c("tread 'solo'", paste0(trees[i], ";"), "proc-;"), file.path(wd, "solo.tre"))
  v <- sapply(1:3, function(s) runbest(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
        "hold 1000;","proc solo.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",30))))
  solo[i] <- min(v)
  cat(sprintf("  tree %2d solo strict best (over seeds): %g   {%s}\n", i, min(v), paste(v,collapse=",")))
}
cat(sprintf("\n  BEST any-single-tree-solo = %g   |  SET-strict reaches 1261/target\n", min(solo)))
cat(sprintf("  => %s\n", if (min(solo) > 1261)
   "no single member reaches target solo: cross-set variety IS the mechanism"
   else "a single member reaches target solo: 'lucky tree', re-examine"))
