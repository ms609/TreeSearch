source("split-support/config.R")

if(!dir.exists("split-support/alignments")) {
  source("split-support/simulate.R")
}

template <- readLines("split-support/tnt-template.run")
on.exit(unlink("split-support/tnt.run"))

for (i in aln) {
  writeLines(gsub("aln####", aln, template, fixed = TRUE),
            "split-support/tnt.run")
  system2(paste0(tntExec, " proc ",  "split-support/tnt.run"))
}
