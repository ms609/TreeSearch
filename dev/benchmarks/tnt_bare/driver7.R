source("dev/benchmarks/tnt_bare/harness.R")
# Independent-parallel ("10 random starts, take best") vs SHARED-buffer (the real set run).
# Same seed (1), same rounds (30), strict noequals throughout. If set < min(independent solos),
# the 10 tracks are NOT independent -- they combine through TNT's single shared tree buffer.
L <- readLines(file.path(bare, paste0(nm, ".t0.tre")))
trees <- sub("[*;]$", "", grep("^\\(", L, value = TRUE))
best <- function(lines) min(rss_bests(run_tnt(c(lines, "quit;"))))
strict <- function(src, h = 1000) c("mxram 1024;","proc data.tnt;","rseed 1;",
           sprintf("hold %d;", h), sprintf("proc %s;", src),
           "sectsch: noglobal noequals;", rep("sectsch=rss;", 30))

# (1) each of the 10 trees SOLO, seed 1, 30 rounds -> independent-parallel baseline
solo <- numeric(length(trees))
for (i in seq_along(trees)) {
  writeLines(c("tread 'solo'", paste0(trees[i], ";"), "proc-;"), file.path(wd, "solo.tre"))
  solo[i] <- best(strict("solo.tre"))
}
cat(sprintf("INDEPENDENT (10 solos, seed1, 30 rnds): each = {%s}\n", paste(solo, collapse=",")))
cat(sprintf("   -> best-of-10 independent = %g\n", min(solo)))

# (2) the 10-tree SET in one shared buffer, seed 1, 30 rounds
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)
set_best <- best(strict("set.tre"))
cat(sprintf("SHARED   (10-tree set, seed1, 30 rnds)  = %g\n", set_best))

cat(sprintf("\nVERDICT: %s (independent-best %g vs shared %g)\n",
    if (set_best < min(solo)) "SHARED beats best-independent -> tracks COMBINE, not just parallel"
    else "shared == best-independent -> consistent with mere parallel restarts",
    min(solo), set_best))
