source("dev/benchmarks/tnt_bare/harness.R")
# Parse helper: pull all RSS best scores + final TreeLength from a raw TNT run that
# writes finalt.tre.
rss_bests <- function(out) num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                  grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
score_final <- function() {            # MIN TreeLength over ALL saved trees (best in memory)
  ff <- file.path(wd, "finalt.tre"); if (!file.exists(ff)) return(NA)
  tr <- tryCatch(ReadTntTree(ff), error = function(e) NULL); if (is.null(tr)) return(NA)
  if (!inherits(tr, "multiPhylo")) tr <- structure(list(tr), class = "multiPhylo")
  tryCatch(min(vapply(tr, function(x) TreeLength(x, phy), numeric(1))), error = function(e) NA)
}
runblk <- function(lines) { out <- run_tnt(c(lines, "tsave *finalt.tre;","save;","tsave/;","quit;"))
  list(bests = rss_bests(out), TL = score_final()) }
file.copy(file.path(bare, paste0(nm, ".t0.tre")), file.path(wd, "set.tre"), overwrite = TRUE)  # 10-tree 1271 set

cat("==== A: SINGLE 1271 tree, various knobs (hold 1000) ====\n")
for (cfg in c("", "noglobal noequals", "equals", "global 1", "equals global 1")) {
  r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","proc tee.tre;",
                if (nzchar(cfg)) sprintf("sectsch: %s;", cfg) else character(0),
                rep("sectsch=rss;", 12)))
  cat(sprintf("  [%-22s] rounds: %s | final TL=%s\n", cfg,
              paste(r$bests, collapse=" "), format(r$TL)))
}

cat("\n==== B: 10-tree 1271 SET start (hold 1000) ====\n")
for (cfg in c("", "noglobal noequals", "equals")) {
  r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1000;","proc set.tre;",
                if (nzchar(cfg)) sprintf("sectsch: %s;", cfg) else character(0),
                rep("sectsch=rss;", 12)))
  cat(sprintf("  [%-22s] rounds: %s | final TL=%s\n", cfg,
              paste(r$bests, collapse=" "), format(r$TL)))
}

cat("\n==== C: in-memory hold-1 mult (=1275) then sectsch (reproduce prior seq_accum) ====\n")
for (cfg in c("", "noglobal noequals", "equals")) {
  r <- runblk(c("mxram 1024;","proc data.tnt;","rseed 1;","hold 1;","mult=replic 1;",
                if (nzchar(cfg)) sprintf("sectsch: %s;", cfg) else character(0),
                rep("sectsch=rss;", 12)))
  cat(sprintf("  [%-22s] rounds: %s | final TL=%s\n", cfg,
              paste(r$bests, collapse=" "), format(r$TL)))
}
