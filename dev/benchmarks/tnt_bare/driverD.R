source("dev/benchmarks/tnt_bare/harness.R")
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
best <- function(lines) min(rss_bests(run_tnt(c(lines, "quit;"))))
hit <- function(v) sprintf("1261 in %d/%d, med=%g [%g-%g]", sum(v<=1261), length(v), median(v), min(v), max(v))

# (1) Does TNT's sector SIZE escalate across sectsch=rss commands? Dump settings between rounds.
cat("==== TNT size schedule: settings between successive sectsch=rss rounds (single tree) ====\n")
out <- run_tnt(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","proc tee.tre;",
                 "sect:;", "sectsch=rss;", "sect:;", "sectsch=rss;", "sect:;", "sectsch=rss;",
                 "sect:;", "quit;"))
size_lines <- grep("size|sectors of|selections", out, ignore.case = TRUE, value = TRUE)
cat(paste0("  ", trimws(size_lines)), sep = "\n")

# (2) Does the escalating schedule matter, or is fixed n/2 enough? SET, fixed vs default size.
cat("\n==== SET (10 trees): fixed sector size vs default(escalating) schedule, seeds 1-6 ====\n")
SEEDS <- 1:6
for (cfg in c("minsize 37 maxsize 37", "minsize 37 maxsize 37 increase 0", "")) {
  v <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
        "hold 1000;","proc set.tre;",
        paste0("sectsch: noglobal noequals", if (nzchar(cfg)) paste0(" ", cfg) else ""), ";",
        rep("sectsch=rss;",30))))
  cat(sprintf("  [%-32s] %s\n", if (nzchar(cfg)) cfg else "default(escalating)", hit(v)))
}
