# Calculate site concordance factor

The site concordance factor (Minh et al. 2020) is a measure of the
strength of support that the dataset presents for a given split in a
tree.

## Usage

``` r
QuartetConcordance(tree, dataset = NULL, weight = TRUE)

ClusteringConcordance(tree, dataset)

PhylogeneticConcordance(tree, dataset)

MutualClusteringConcordance(tree, dataset)

SharedPhylogeneticConcordance(tree, dataset)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree. Perhaps load into R
  using
  [`ReadAsPhyDat()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html).
  Additive (ordered) characters can be handled using
  [`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html).

- weight:

  Logical specifying whether to weight sites according to the number of
  quartets they are decisive for.

## Details

`QuartetConcordance()` is the proportion of quartets (sets of four
leaves) that are decisive for a split which are also concordant with it.
For example, a quartet with the characters `0 0 0 1` is not decisive, as
all relationships between those leaves are equally parsimonious. But a
quartet with characters `0 0 1 1` is decisive, and is concordant with
any tree that groups the first two leaves together to the exclusion of
the second.

By default, the reported value weights each site by the number of
quartets it is decisive for. This value can be interpreted as the
proportion of all decisive quartets that are concordant with a split. If
`weight = FALSE`, the reported value is the mean of the concordance
value for each site. Consider a split associated with two sites: one
that is concordant with 25% of 96 decisive quartets, and a second that
is concordant with 75% of 4 decisive quartets. If `weight = TRUE`, the
split concordance will be 24 + 3 / 96 + 4 = 27%. If `weight = FALSE`,
the split concordance will be mean(75%, 25%) = 50%.

`QuartetConcordance()` is computed exactly, using all quartets, where as
other implementations (e.g. IQ-TREE) follow Minh2020) in using a random
subsample of quartets for a faster, if potentially less accurate,
computation.

**NOTE:** These functions are under development. They are incompletely
tested, and may change without notice. Complete documentation and
discussion will follow in due course.

## References

Minh BQ, Hahn MW, Lanfear R (2020). “New methods to calculate
concordance factors for phylogenomic datasets.” *Molecular Biology and
Evolution*, **37**(9), 2727–2733.
[doi:10.1093/molbev/msaa106](https://doi.org/10.1093/molbev/msaa106) .

## See also

Other split support functions:
[`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md),
[`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md),
[`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md),
[`MostContradictedFreq()`](https://ms609.github.io/TreeSearch/reference/MostContradictedFreq.md),
[`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("congreveLamsdellMatrices", package = "TreeSearch")
dataset <- congreveLamsdellMatrices[[1]][, 1:20]
tree <- referenceTree
qc <- QuartetConcordance(tree, dataset)
cc <- ClusteringConcordance(tree, dataset)
pc <- PhylogeneticConcordance(tree, dataset)
spc <- SharedPhylogeneticConcordance(tree, dataset)
mcc <- MutualClusteringConcordance(tree, dataset)

oPar <- par(mar = rep(0, 4), cex = 0.8) # Set plotting parameters
plot(tree)
TreeTools::LabelSplits(tree, signif(qc, 3), cex = 0.8)

plot(tree)
TreeTools::LabelSplits(tree, signif(cc, 3), cex = 0.8)

par(oPar) # Restore plotting parameters

# Write concordance factors to file
labels <- paste0(qc, "/", cc, "/", pc) # "/" is a valid delimiter
# Identify the node that corresponds to each label
whichNode <- match(TreeTools::NTip(tree) + 1:tree$Nnode, names(qc))

# The contents of tree$node.label will be written at each node
tree$node.label <- labels[whichNode]

ape::write.tree(tree) # or write.nexus(tree, file = "mytree.nex")
#> [1] "(((((((((((((2,3)0.794157608695652/0.0840869199349208/0.8,4)0.796019900497512/0.106844441779741/0.7,5)0.745984310795667/0.125579285095439/0.5,6)0.744573765335011/0.157227405925398/0.45,7)0.726778007130794/0.158079354442316/0.45,8)0.72615039281706/0.161803744635535/0.45,9)0.736253135249855/0.164460425292521/0.45,10)0.713525905547914/0.152606841335018/0.35,11)0.693087173792338/0.139772020193245/0.35,12)0.663097646420289/0.119374128613331/0.35,(19,((18,(17,(15,16)0.575/0.0480987344280996/0.75)0.62020202020202/0.0630658263101711/0.65)0.606779661016949/0.0675277845468592/0.5,(13,14)0.6825/0.0463478996057617/0.85)0.628180398640581/0.0893831815805849/0.45)0.649913605288859/0.0972690748771144/0.4)0.665/0.0733038953664147/0.5,(20,(21,22)0.631908237747654/0.0498507296733025/0.8)0.60419091967404/0.0567591383942297/0.65)NA,1)NA;"

# Display correlation between concordance factors
pairs(cbind(qc, cc, pc, spc, mcc), asp = 1)
```
