# When was a tree topology first hit?

Reports when each tree in a list was first found by tree search. This
information is read from the `firstHit` attribute if present. If not,
trees are taken to be listed in the order in which they were found, and
named according to the search iteration in which they were first hit -
the situation when trees found by
[`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)
are saved to file.

## Usage

``` r
WhenFirstHit(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object.

## Value

`trees`, with a `firstHit` attribute listing the number of trees hit for
the first time in each search iteration.

## See also

- [`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md)

Other utility functions:
[`ClusterStrings()`](https://ms609.github.io/TreeSearch/reference/ClusterStrings.md),
[`QACol()`](https://ms609.github.io/TreeSearch/reference/QACol.md),
[`QuartetResolution()`](https://ms609.github.io/TreeSearch/reference/QuartetResolution.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("TreeTools", quietly = TRUE)
trees <- list(
   seed_00 = as.phylo(1, 8),
   ratch1_01 = as.phylo(2, 8),
   ratch1_02 = as.phylo(3, 8),
   ratch4_44 = as.phylo(4, 8),
   final_99 = as.phylo(5, 8)
)
attr(WhenFirstHit(trees), "firstHit")
#> whenHit
#>   seed ratch1 ratch4  final 
#>      1      2      1      1 
```
