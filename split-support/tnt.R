if(!dir.exists("split-support/alignments")) {
  source("simulate.R")
}

# Download command-driven TNT from 
# http://www.lillo.org.ar/phylogeny/tnt/ZIPCHTNT.ZIP
# and unzip to a convenient local folder.

tntPath <- "C:/Programs/Phylogeny/tnt/tnt.exe"

template <- readLines("split-support/tnt-template.run")
on.exit(unlink("split-support/tnt.run"))

for (i in 1:1000) {
  writeLines(gsub("####", formatC(i, width = 4, flag = 0), template),
            "split-support/tnt.run")
  system2(paste0(tntPath, " proc ",  "split-support/tnt.run"))
}

