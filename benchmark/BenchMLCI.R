library("TreeTools", quietly = TRUE)
library("TreeSearch")

set.seed(0)
tree <- referenceTree
data <- PhyDatToMatrix(congreveLamsdellMatrices[[1]])

smallTree <- KeepTip(referenceTree, as.character(1:9))
smallData <- data[as.character(1:9), ]

bench::mark(
  small = MLCI(smallTree, smallData, precision = 1/20),
  medium = MLCI(tree, data, precision = 1/20),
  max_iterations = 12,
  check = FALSE
)

bench::mark(
  small = MLCI(smallTree, smallData, precision = 1/20),
  medium = MLCI(tree, data, precision = 1/20),
  min_time = 22,
  check = FALSE
)
