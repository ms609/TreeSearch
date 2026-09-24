source("dev/benchmarks/tnt_bare/harness.R")
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
runblk <- function(lines) { out <- run_tnt(c(lines, "tsave *finalt.tre;","save;","tsave/;","quit;"))
  list(bests = rss_bests(out), TL = score_final(), n = n_trees()) }

cat("==== CHECK 1: TNT's ACTUAL DEFAULT pipeline (noequals throughout) ====\n")
cat("     proc; rseed 1; hold H; mult=replic 1;  then sectsch=rss to plateau\n")
for (h in c(1, 1000)) {
  r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;",sprintf("hold %d;",h),
                "mult=replic 1;", rep("sectsch=rss;",16)))
  cat(sprintf("  hold=%-4d : best=%s  ntrees=%d  rounds: %s\n", h, format(min(r$bests)),
              r$n, paste(r$bests, collapse=" ")))
}

cat("\n==== CHECK 2: variety accumulation, single T0 (A): trees-in-memory after k rounds ====\n")
for (cfg in c("noglobal noequals", "equals")) {
  cat(sprintf("  -- %s --\n", cfg))
  for (k in c(1,2,4,8,12)) {
    r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","proc tee.tre;",
                  sprintf("sectsch: %s;",cfg), rep("sectsch=rss;",k)))
    cat(sprintf("     k=%2d : best=%s  ntrees=%d\n", k, format(min(r$bests)), r$n))
  }
}

cat("\n==== CHECK 3: SET + STRICT — is it fusing or multi-tree sectorial? ====\n")
for (cfg in c("noglobal noequals", "noglobal noequals nofuse", "noglobal noequals tree 0")) {
  r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","proc set.tre;",
                sprintf("sectsch: %s;",cfg), rep("sectsch=rss;",12)))
  cat(sprintf("  [%-26s] best=%s  ntrees=%d  rounds: %s\n", cfg, format(min(r$bests)),
              r$n, paste(r$bests, collapse=" ")))
}
