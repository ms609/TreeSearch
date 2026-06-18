# BARE-BONES sectsch from a fixed single-tree T0, read fresh (NO mult/TBR before sectsch).
# Dumps raw TNT output so we can see every reported score / accepted move.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"), winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
bare <- "dev/benchmarks/tnt_bare"
nm <- Sys.getenv("DS", "Zanol2014")
phy <- readRDS(file.path(bare, paste0(nm, ".phy.rds")))

wd <- file.path(tempdir(), paste0("bb", Sys.getpid(), nm)); unlink(wd, recursive = TRUE)
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))

# ---- Step 1: build a SINGLE-tree T0 (mult replic 1, rseed 1, hold 1), save to t0.tre ----
writeLines(c("mxram 1024;", "proc data.tnt;", "rseed 1;", "hold 1;",
             "mult = replic 1;", "tsave *tee.tre;", "save;", "tsave/;", "quit;"),
           file.path(wd, "maketee.run"))
old <- setwd(wd)
suppressWarnings(system2(TNT, args = "maketee.run;", stdout = TRUE, stderr = TRUE))
setwd(old)
t0 <- ReadTntTree(file.path(wd, "tee.tre")); if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
cat(sprintf("T0 single-tree score = %.0f (tips=%d)\n", TreeLength(t0, phy), length(t0$tip.label)))
file.copy(file.path(wd, "tee.tre"), file.path(bare, paste0(nm, ".t0single.tre")), overwrite = TRUE)

# ---- Step 2: FRESH session: load matrix, read T0, run ONE stripped sectsch round ----
script <- Sys.getenv("SCRIPT_FILE", "")
if (!nzchar(script)) {
  script <- file.path(wd, "barerun.run")
  writeLines(c(
    "mxram 1024;",
    "report+;",                 # verbose progress
    "proc data.tnt;",
    "rseed 1;",
    "hold 1000;",
    "proc tee.tre;",            # read the fixed T0 (no search yet)
    "sect: ;",                  # show CURRENT (default) sectsch settings
    # strip the obvious bells: no global TBR, strict acceptance, no fusing, no drift
    "sectsch: noglobal noequals nofuse godrift 9999 ;",
    "sect: ;",                  # show stripped settings
    "sectsch = rss ;",          # ONE bare sectorial round
    "score ;",
    "quit;"), script)
} else {
  file.copy(script, file.path(wd, basename(script)), overwrite = TRUE)
  script <- file.path(wd, basename(script))
}
old <- setwd(wd)
out <- suppressWarnings(system2(TNT, args = paste0(basename(script), ";"), stdout = TRUE, stderr = TRUE))
setwd(old)
out <- iconv(out, from = "", to = "UTF-8", sub = "")
cat("==== RAW TNT OUTPUT (filtered to informative lines) ====\n")
keep <- grep("score|RSS|ector|eplac|earrang|settings|size|global|equal|drift|fuse|RAS|TBR",
             out, ignore.case = TRUE, value = TRUE)
cat(paste0(trimws(keep)), sep = "\n")
