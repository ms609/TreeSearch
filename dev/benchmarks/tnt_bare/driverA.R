source("dev/benchmarks/tnt_bare/harness.R")
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
trees <- sub("[*;]$", "", grep("^\\(", L, value = TRUE))
best <- function(lines) min(rss_bests(run_tnt(c(lines, "quit;"))))
strict_rounds <- function(src, s, R = 30) c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
   "hold 1000;", sprintf("proc %s;", src), "sectsch: noglobal noequals;", rep("sectsch=rss;", R))

# ---- TEST 1: can ABUNDANT independent restarts reach 1261? (is it 'just restarts'?) ----
# 10 trees x 20 seeds = 200 independent single-tree strict runs.
NS <- 20
allbest <- 9999; hit <- 0; n <- 0
for (i in seq_along(trees)) {
  writeLines(c("tread 'solo'", paste0(trees[i], ";"), "proc-;"), file.path(wd, "solo.tre"))
  for (s in 1:NS) { b <- best(strict_rounds("solo.tre", s)); allbest <- min(allbest, b)
    hit <- hit + (b <= 1261); n <- n + 1 }
}
cat(sprintf("TEST1  %d independent solo restarts: best=%g, #reaching<=1261 = %d/%d\n", n, allbest, hit, n))

# ---- TEST 2: isolate COUPLING from diversity. Same single tree, 10 copies. ----
# (a) 10 copies in ONE shared buffer (strict) vs (b) 10 independent restarts of that tree.
ti <- 1  # tree1 is frozen solo at seed1
writeLines(c("tread 'solo'", paste0(trees[ti], ";"), "proc-;"), file.path(wd, "solo.tre"))
ident <- c(L[1], paste0(rep(paste0(trees[ti], "*"), 9), collapse = "\n"), paste0(trees[ti], ";"), "proc-;")
writeLines(ident, file.path(wd, "ident.tre"))
for (s in 1:5) {
  shared <- best(strict_rounds("ident.tre", s))        # 10 identical copies, shared buffer
  indep  <- min(sapply(1:10, function(ss) best(strict_rounds("solo.tre", (s-1)*10+ss))))  # 10 separate restarts
  cat(sprintf("TEST2  seed-block %d : 10-copies-shared-buffer=%g   vs  10-separate-restarts=%g\n",
              s, shared, indep))
}
