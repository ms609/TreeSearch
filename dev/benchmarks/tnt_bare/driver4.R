source("dev/benchmarks/tnt_bare/harness.R")
runblk <- function(lines) { out <- run_tnt(c(lines, "tsave *finalt.tre;","save;","tsave/;","quit;"))
  list(bests = rss_bests(out), TL = score_final(), n = n_trees(),
       err = grep("rror|nvalid|nrecogni", out, value = TRUE)) }

# ---- Build a 10-IDENTICAL-copies set (tree 0 x10) to control for tree COUNT vs DIVERSITY ----
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
tree0 <- sub("[*]$", "", L[2])
ident <- c(L[1], paste0(rep(paste0(tree0, "*"), 9), collapse = "\n"), paste0(tree0, ";"), "proc-;")
writeLines(ident, file.path(wd, "ident.tre"))
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)

cat("==== VARIETY CONTROL: 10 IDENTICAL copies vs 10 DIFFERENT trees (strict sectsch) ====\n")
for (src in c("ident.tre", "set.tre")) {
  r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;",sprintf("proc %s;",src),
                "sectsch: noglobal noequals;", rep("sectsch=rss;",12)))
  cat(sprintf("  [%-10s] best=%s ntrees(end)=%d  rounds: %s\n", src, format(min(r$bests)),
              r$n, paste(r$bests, collapse=" ")))
}

cat("\n==== `tree 0` behaviour check (did it error / restrict?) ====\n")
r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","proc set.tre;",
              "sectsch: noglobal noequals tree 0;", rep("sectsch=rss;",3)))
cat(sprintf("  errors: %s\n", if (length(r$err)) paste(unique(trimws(r$err)),collapse=" | ") else "<none>"))

cat("\n==== SEED ROBUSTNESS ====\n")
cat("-- Canonical default pipeline (mult + sectsch, noequals), hold 1 vs 1000 --\n")
for (h in c(1, 1000)) for (s in 1:3) {
  r <- runblk(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),sprintf("hold %d;",h),
                "mult=replic 1;", rep("sectsch=rss;",16)))
  cat(sprintf("  hold=%-4d seed=%d : best=%s ntrees=%d\n", h, s, format(min(r$bests)), r$n))
}
cat("-- Single T0 (fixed seed-1 tree), vary sectsch rseed --\n")
for (cfg in c("noglobal noequals", "equals")) for (s in 1:3) {
  r <- runblk(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),"hold 1000;","proc tee.tre;",
                sprintf("sectsch: %s;",cfg), rep("sectsch=rss;",12)))
  cat(sprintf("  [%-18s] sectsch-seed=%d : best=%s ntrees=%d\n", cfg, s, format(min(r$bests)), r$n))
}
