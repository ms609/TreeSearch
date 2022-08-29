source("split-support/config.R")

if(!dir.exists("split-support/alignments")) {
  source("split-support/simulate.R")
}

template <- readLines("split-support/tnt-template.run")
on.exit(unlink("split-support/tnt.run"))

for (aln in alns) {
  if (file.exists(paste0("split-support/tnt/", aln, ".iw16.sym"))) {
    message("Results found for ", aln)
  } else {
    writeLines(gsub("aln####", aln, template, fixed = TRUE),
              "split-support/tnt.run")
    system2(tntExec, "proc split-support/tnt.run")
  }
}

warnings()
