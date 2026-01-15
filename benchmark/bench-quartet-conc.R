library("TreeTools")

tmp_lib <- tempfile(pattern = "lib")
dir.create(tmp_lib)
devtools::install(args = paste0("--library=", tmp_lib))
library("TreeSearch", lib.loc = tmp_lib)

data("congreveLamsdellMatrices", package = "TreeSearch")
dataset <- congreveLamsdellMatrices[[42]]
someNA <- PhyDatToMatrix(dataset)
someNA[sample.int(length(someNA), 444)] <- "?"
someNA <- MatrixToPhyDat(someNA)
tree <- TreeSearch::referenceTree

bench::mark(
  QuartetConcordance(tree, dataset),
  QuartetConcordance(tree, someNA), check = FALSE
)
