library("TreeTools")
library("TreeSearch")
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

#
# A tibble: 2 × 13
  # expression         min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory    
  # <bch:expr>       <bch> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list>    
# 1 QuartetConcorda… 237ms  239ms      4.19    2.01MB     2.09     2     1      478ms <NULL> <Rprofmem>
# 2 QuartetConcorda… 233ms  236ms      4.23    1.88MB     2.12     2     1      472ms <NULL> <Rprofmem>
# ℹ 2 more variables: time <list>, gc <list>
