# Dump TNT help for the commands relevant to Wagner-only (no-swap) RAS and
# fixed-addition-sequence, plus bbreak/mult, to verify exact syntax before
# spending the K=200 TNT batch.  Uses the define_target.R system2 pattern.
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
wd <- file.path(tempdir(), paste0("tnthelp", Sys.getpid()))
unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  "mxram 1024;",
  "help mult;",
  "help rseed;",
  "help bbreak;",
  "help randtrees;",
  "help hold;",
  "quit;"),
  file.path(wd, "helpdump.run"))
old <- setwd(wd)
out <- suppressWarnings(system2(TNT, args = "helpdump.run;", stdout = TRUE, stderr = TRUE))
setwd(old)
out <- iconv(out, from = "", to = "UTF-8", sub = "")
cat(out, sep = "\n")
