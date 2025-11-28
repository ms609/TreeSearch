library("TreeTools", quietly = TRUE)
library("TreeSearch")

set.seed(0)
tree <- referenceTree
data <- PhyDatToMatrix(congreveLamsdellMatrices[[1]])

smallTree <- KeepTip(referenceTree, as.character(1:9))
smallData <- data[as.character(1:9), ]

bench::mark(
  small = MLCI(smallTree, smallData, precision = 1/10),
  medium = MLCI(tree, data, precision = 1/10),
  max_iterations = 12,
  check = FALSE
)

bench::mark(
  single = {MaddisonSlatkin_clear_cache(); MaddisonSlatkin(3L, states)},
  check = FALSE,
  min_time = 20
)
