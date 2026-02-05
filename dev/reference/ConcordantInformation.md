# Evaluate the concordance of information between a tree and a dataset

Details the amount of information in a phylogenetic dataset that is
consistent with a specified phylogenetic tree, and the signal:noise
ratio of the character matrix implied if the tree is true.

## Usage

``` r
ConcordantInformation(tree, dataset)

Evaluate(tree, dataset)

ConcordantInfo(tree, dataset)
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

## Value

`ConcordantInformation()` returns a named vector with elements:

- `informationContent`: cladistic information content of `dataset`

- `signal`, `noise`: amount of cladistic information that represents
  phylogenetic signal and noise, according to `tree`

- `signalToNoise`: the implied signal:noise ratio of `dataset`

- `treeInformation`: the cladistic information content of a bifurcating
  tree on `dataset`; this is the minimum amount of information necessary
  to resolve a bifurcating tree, assuming no duplicate information or
  noise

- `matrixToTree`: the ratio of the cladistic information content of the
  matrix to the cladistic information content of the tree, a measure of
  the redundancy of the matrix

- `ignored`: information content of characters whose signal and noise
  could not be calculated (too many states) and so are not included in
  the totals above.

## Details

Presently restricted to datasets whose characters contain a maximum of
two parsimony-informative states.

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data(congreveLamsdellMatrices)
myMatrix <- congreveLamsdellMatrices[[10]]
ConcordantInformation(TreeTools::NJTree(myMatrix), myMatrix)
#> dataset contains 821.038 bits, of which 280.813 signal, 540.226 noise, 78.0817 needed.  S:N = 0.519806
#> informationContent             signal              noise      signalToNoise 
#>        821.0383923        280.8125010        540.2258914          0.5198057 
#>    treeInformation       matrixToTree            ignored 
#>         78.0816559         10.5151253          0.0000000 
```
