# Addition tree

Generates a starting tree by adding each taxon in turn to the most
parsimonious location.

## Usage

``` r
AdditionTree(dataset, concavity = Inf, constraint, sequence)
```

## Arguments

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree. Perhaps load into R
  using
  [`ReadAsPhyDat()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html).
  Additive (ordered) characters can be handled using
  [`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html).

- concavity:

  Determines the degree to which extra steps beyond the first are
  penalized. Specify a numeric value to use implied weighting
  (Goloboff 1993) ; `concavity` specifies *k* in *k* / *e* + *k*. A
  value of 10 is recommended; TNT sets a default of 3, but this is too
  low in some circumstances (Goloboff et al. 2018; Smith 2019) . Better
  still explore the sensitivity of results under a range of concavity
  values, e.g. `k = 2 ^ (1:7)`. Specify `Inf` to weight each additional
  step equally, (which underperforms step weighting approaches (Goloboff
  et al. 2008; Goloboff et al. 2018; Goloboff and Arias 2019;
  Smith 2019) ). Specify `"profile"` to employ an approximation of
  profile parsimony (Faith and Trueman 2001) .

- constraint:

  Either an object of class `phyDat`, in which case returned trees will
  be perfectly compatible with each character in `constraint`; or a tree
  of class `phylo`, all of whose nodes will occur in any output tree.
  See
  [`ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint.html)
  and
  [vignette](https://ms609.github.io/TreeSearch/articles/tree-search.html)
  for further examples.

- sequence:

  Character or numeric vector listing sequence in which to add taxa.
  Randomized if not provided.

## Value

`AdditionTree()` returns a tree of class `phylo`, rooted on
`sequence[1]`.

## See also

Impose a constraint:
[[`TreeTools::ImposeConstraint()`](https://ms609.github.io/TreeTools/reference/ImposeConstraint.html)](https://ms609.github.io/TreeTools/reference/ImposeConstraint)

Neighbour-joining trees:
[[`TreeTools::NJTree()`](https://ms609.github.io/TreeTools/reference/NJTree.html)](https://ms609.github.io/TreeTools/reference/NJTree.html);
[[`TreeTools::ConstrainedNJ()`](https://ms609.github.io/TreeTools/reference/ConstrainedNJ.html)](https://ms609.github.io/TreeTools/reference/ConstrainedNJ)

Other tree generation functions:
[`RandomMorphyTree()`](https://ms609.github.io/TreeSearch/reference/RandomMorphyTree.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("inapplicable.phyData", package = "TreeSearch")
AdditionTree(inapplicable.phyData[["Longrich2010"]], concavity = 10)
#> 
#> Phylogenetic tree with 20 tips and 19 internal nodes.
#> 
#> Tip labels:
#>   Psittacosaurus_spp, Pachycephalosaurus_wyomingensis, Texacephale_langstoni, Hanssuesia_sternbergi, Homalocephale_calathocercos, Goyocephale_lattimorei, ...
#> 
#> Rooted; no branch length.
```
