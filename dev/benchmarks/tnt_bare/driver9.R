source("dev/benchmarks/tnt_bare/harness.R")
# Is the set just "restart from 10 diverse trees, take best" (no sharing), or does the shared
# population reach something independent restarts cannot? Give each tree its OWN seeds (proper
# independent-restart baseline) and pour on restarts; compare best to the set's 1261.
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
trees <- sub("[*;]$", "", grep("^\\(", L, value = TRUE))
best <- function(lines) min(rss_bests(run_tnt(c(lines, "quit;"))))
solo_best <- function(i, s) { writeLines(c("tread 'solo'", paste0(trees[i], ";"), "proc-;"),
    file.path(wd, "solo.tre"))
  best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),"hold 1000;","proc solo.tre;",
         "sectsch: noglobal noequals;", rep("sectsch=rss;",30))) }

SEEDS <- 1:5
M <- outer(seq_along(trees), SEEDS, Vectorize(function(i,s) solo_best(i,s)))
rownames(M) <- paste0("tree", seq_along(trees)); colnames(M) <- paste0("s", SEEDS)
cat("Per-tree solo strict best (rows=tree, cols=seed):\n")
print(M)
diag_seed <- sapply(seq_along(trees), function(i) M[i, ((i-1) %% length(SEEDS))+1])
cat(sprintf("\n'10 independent restarts (tree i, its own seed)' best = %g\n", min(diag_seed)))
cat(sprintf("BEST over ALL %d independent solo runs                = %g\n", length(M), min(M)))
cat(sprintf("Number of independent runs that reached <=1261        = %d / %d\n",
            sum(M <= 1261), length(M)))

# the shared set, matched seeds
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
set_best <- sapply(SEEDS, function(s) best(c("mxram 1024;","proc data.tnt;",sprintf("rseed %d;",s),
   "hold 1000;","proc set.tre;","sectsch: noglobal noequals;", rep("sectsch=rss;",30))))
cat(sprintf("\nSHARED 10-tree set best over seeds 1-5                 = %g  {%s}\n",
            min(set_best), paste(set_best, collapse=",")))
cat(sprintf("\nVERDICT: %s\n", if (min(M) <= min(set_best))
  "independent restarts MATCH the set -> it's 'restart from diverse trees, take best', NOT sharing"
  else "set BEATS all independent restarts -> genuine population synergy / sharing"))
