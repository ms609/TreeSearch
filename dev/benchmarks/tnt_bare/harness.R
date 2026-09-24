# Reusable TNT sectsch harness.
#   Reads a FIXED single-tree T0 fresh (no mult/TBR before sectsch),
#   applies a sectsch config, runs N rounds of `sectsch=rss;`, capturing the
#   running best score after each round + the final score. TNT score is
#   authoritative; the final tree is also re-scored with TreeLength as a
#   mapping sanity check.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"), winslash = "/"))
  library(TreeTools)
})
TNT  <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
bare <- "dev/benchmarks/tnt_bare"
nm   <- Sys.getenv("DS", "Zanol2014")
phy  <- readRDS(file.path(bare, paste0(nm, ".phy.rds")))
t0file <- file.path(bare, paste0(nm, ".t0single.tre"))
num  <- function(x) suppressWarnings(as.double(gsub(",", "", x)))

# One reusable working dir per process
wd <- file.path(tempdir(), paste0("hn", Sys.getpid(), nm)); unlink(wd, recursive = TRUE)
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))
file.copy(t0file, file.path(wd, "tee.tre"), overwrite = TRUE)

rss_bests <- function(out) num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                  grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
read_trees <- function(ff) {
  if (!file.exists(ff)) return(NULL)
  tr <- tryCatch(ReadTntTree(ff), error = function(e) NULL); if (is.null(tr)) return(NULL)
  if (!inherits(tr, "multiPhylo")) tr <- structure(list(tr), class = "multiPhylo"); tr
}
score_final <- function(ff = file.path(wd, "finalt.tre")) {  # MIN TreeLength over saved trees
  tr <- read_trees(ff); if (is.null(tr)) return(NA)
  tryCatch(min(vapply(tr, function(x) TreeLength(x, phy), numeric(1))), error = function(e) NA)
}
n_trees <- function(ff = file.path(wd, "finalt.tre")) { tr <- read_trees(ff); if (is.null(tr)) NA else length(tr) }

run_tnt <- function(lines) {
  rf <- file.path(wd, "runme.run")
  writeLines(lines, rf)
  old <- setwd(wd)
  out <- suppressWarnings(system2(TNT, args = "runme.run;", stdout = TRUE, stderr = TRUE))
  setwd(old)
  iconv(out, from = "", to = "UTF-8", sub = "")
}

# Run a config: setting_line is the `sectsch: ...;` options (may be ""), rounds = #sectsch=rss
run_config <- function(setting_line, rounds = 8, seed = 1, hold = 1000, label = "") {
  pre <- c("mxram 1024;", "proc data.tnt;", sprintf("rseed %d;", seed),
           sprintf("hold %d;", hold), "proc tee.tre;")
  if (nzchar(setting_line)) pre <- c(pre, sprintf("sectsch: %s;", setting_line))
  body <- as.vector(rbind(rep("sectsch = rss;", rounds),
                          rep("tplot/;", 0)))           # placeholder, removed below
  body <- rep("sectsch = rss;", rounds)
  lines <- c(pre, "score;", body, "score;", "tsave *finalt.tre;", "save;", "tsave/;", "quit;")
  out <- run_tnt(lines)
  # Parse every "best score:" from RSS, plus start/end "score;" outputs
  best_lines <- grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)
  bests <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1", best_lines))
  # final tree score via TreeLength (mapping check)
  tl <- NA
  ff <- file.path(wd, "finalt.tre")
  if (file.exists(ff)) {
    tr <- tryCatch(ReadTntTree(ff), error = function(e) NULL)
    if (!is.null(tr)) { if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
      tl <- tryCatch(min(TreeLength(tr, phy)), error = function(e) NA) }
  }
  list(label = label, setting = setting_line, rounds = rounds, seed = seed, hold = hold,
       per_round = bests, final_tnt = if (length(bests)) min(bests) else NA, final_TL = tl)
}

print_config <- function(r) {
  cat(sprintf("\n[%s] hold=%d seed=%d  '%s'\n", r$label, r$hold, r$seed, r$setting))
  cat(sprintf("  per-round best: %s\n", paste(r$per_round, collapse = " ")))
  cat(sprintf("  FINAL TNT=%s  TreeLength=%s\n", format(r$final_tnt), format(r$final_TL)))
}
