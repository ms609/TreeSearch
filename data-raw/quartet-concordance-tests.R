library("TreeTools")
tree <- ape::read.tree(text = "(a, (b, (c, (d, ((e, f), (g, h))))));")
mataset <- matrix(c(0, 0, 0, 0, 1, 1, 1, 1,  0,
                    0, 0, 1, 1, 1, 1, 1, 1,  0,
                    0, 0, 1, 1, 1, 1, 2, 2,  0,
                    1, 0, 0, 0, 1, 1, 1, 1,  0,
                    1, 0, 0, 0, 0, 1, 1, 1,  0,
                    0, 0, 0, 0, 1, 1, 2, 2,  0,
                    0, 0, 1, 1, 2, 2, 3, 3,  0,
                    0, 1, 2, 3, 0, 1, 2, 3,  0,
                    0, 1, 2, 0, 0, 2, 2, 3,  0), 9,
                  dimnames = list(letters[1:9], NULL))
dat <- MatrixToPhyDat(mataset[1:8, ])
iqtreeDir <- stop("Path to folder that contains iqtree2.exe")
treeFile <- paste0(iqtreeDir, "tree.nwk")
dataFile <- paste0(iqtreeDir, "data.nex")
charFile <- paste0(iqtreeDir, "data%i.nex")

write.tree(tree, file = treeFile)

# Expected output per character
vapply(1:9, function(i) {
  charIFile <- sprintf(charFile, i)
  write.nexus.data(dat[, i], file = charIFile, format = "STANDARD")
  dataLines <- readLines(charIFile)
  dataLines[[5]] <- sub("MISSING=? GAP=- INTERLEAVE=NO ", "", dataLines[[5]],
                        fixed = TRUE)
  writeLines(dataLines, charIFile)
  system2(paste0(iqtreeDir, "iqtree2.exe"),
          paste("-t", treeFile, "-s", charIFile,
                "--scf 100000",
                "-nt 4"
          ), stdout = NULL)
  cf <- read.tree(paste0(treeFile, ".cf.tree"))
  
  as.numeric(cf[["node.label"]][-(1:2)])
}, double(5))

# Expected output across entire dataset
write.nexus.data(dat, file = dataFile, format = "STANDARD")
dataLines <- readLines(dataFile)
dataLines[[5]] <- sub("MISSING=? GAP=- INTERLEAVE=NO ", "", dataLines[[5]], fixed = TRUE)
writeLines(dataLines, dataFile)
system2(paste0(iqtreeDir, "iqtree2.exe"),
        paste("-t", treeFile, "-s", dataFile,
              "--scf 1000000",
              "-nt 4"
        ), stdout = NULL)
cf <- read.tree(paste0(treeFile, ".cf.tree"))
plot(cf, show.node.label = TRUE, use.edge.length = FALSE)

as.numeric(cf[["node.label"]][-(1:2)])
