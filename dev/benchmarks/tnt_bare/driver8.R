source("dev/benchmarks/tnt_bare/harness.R")
# Trace the BUFFER across rounds: how many trees, and the spread of their lengths.
# Distinguishes (a) 10 slots independently descending vs (b) buffer collapsing to best &
# re-seeding. Strict noequals, seed 1.
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
writeLines(c("tread 'solo'", paste0(sub("[*;]$","",grep("^\\(",L,value=TRUE))[1],";"), "proc-;"),
           file.path(wd, "solo.tre"))

trace_buffer <- function(src, ks = c(0,1,2,4,8,16,30)) {
  for (k in ks) {
    out <- run_tnt(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;",
                     sprintf("proc %s;", src),
                     "sectsch: noglobal noequals;",
                     if (k > 0) rep("sectsch=rss;", k) else character(0),
                     "tsave *ft.tre;","save;","tsave/;","quit;"))
    tr <- read_trees(file.path(wd, "ft.tre"))
    sc <- if (is.null(tr)) NA else vapply(tr, function(x) TreeLength(x, phy), numeric(1))
    cat(sprintf("  k=%2d : ntrees=%2d  lengths: min=%g max=%g  distinct={%s}\n",
                k, length(sc), min(sc), max(sc), paste(sort(unique(sc)), collapse=",")))
  }
}
cat("==== SET (10 diverse trees) buffer trace ====\n"); trace_buffer("set.tre")
cat("\n==== SOLO (1 tree) buffer trace ====\n");        trace_buffer("solo.tre")
