if(!dir.exists("split-support/alignments")) {
  source("simulate.R")
}

# Download MrBayes from 
# https://nbisweden.github.io/MrBayes/download.html
# and unzip/install to a convenient local folder.

mbPath <-  "C:/Programs/Phylogeny/MrBayes/bin/mb.3.2.7-win64.exe"

template <- readLines("split-support/mb.nex")
on.exit(unlink("split-support/mb.run"))

for (i in formatC(1:1000, width = 4, flag = 0)) {
  writeLines(c(readLines(paste0("split-support/alignments/aln", i, ".nex")),
               template),
             "split-support/mb.run")
  system2(paste0(mbPath, " split-support/mb.run"))
}
