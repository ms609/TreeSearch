source("dev/benchmarks/tnt_bare/harness.R")
# DECISIVE paired test of the user's model: "10 lanes, sectorial within each, pick best."
# For each replicate s: SET (10 trees together, 30 rds) vs 10-INDEPENDENT (same 10 trees, each
# solo at a distinct seed, take min). Equal compute (both = 10 trees x 30 rds). Many replicates.
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
trees <- sub("[*;]$", "", grep("^\\(", L, value = TRUE))
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
best <- function(lines) min(rss_bests(run_tnt(c(lines, "quit;"))))
strict <- function(src, s) c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),"hold 1000;",
   sprintf("proc %s;", src), "sectsch: noglobal noequals;", rep("sectsch=rss;", 30))

REPS <- 15
set_min <- indep_min <- numeric(REPS)
for (s in 1:REPS) {
  set_min[s] <- best(strict("set.tre", s))
  # 10 independent lanes: tree i at its own distinct seed; take the best (min) lane
  solo <- numeric(length(trees))
  for (i in seq_along(trees)) {
    writeLines(c("tread 'solo'", paste0(trees[i], ";"), "proc-;"), file.path(wd, "solo.tre"))
    solo[i] <- best(strict("solo.tre", (s-1)*length(trees) + i))
  }
  indep_min[s] <- min(solo)
  cat(sprintf("  rep %2d : SET=%g   10-INDEP-min=%g   %s\n", s, set_min[s], indep_min[s],
              if (set_min[s] < indep_min[s]) "set<indep" else if (set_min[s] > indep_min[s]) "set>indep" else "tie"))
}
cat(sprintf("\nSET     reached 1261 in %d/%d reps   (median %g)\n", sum(set_min<=1261), REPS, median(set_min)))
cat(sprintf("10-INDEP reached 1261 in %d/%d reps  (median %g)\n", sum(indep_min<=1261), REPS, median(indep_min)))
cat(sprintf("paired: set<indep %d, tie %d, set>indep %d\n",
            sum(set_min<indep_min), sum(set_min==indep_min), sum(set_min>indep_min)))
