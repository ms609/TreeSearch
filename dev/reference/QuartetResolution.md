# Relationship between four taxa

Relationship between four taxa

## Usage

``` r
QuartetResolution(trees, tips)
```

## Arguments

- trees:

  A list of trees of class `phylo`, or a `multiPhylo` object.

- tips:

  Vector specifying four tips whose relationship should be reported, in
  a format accepted by
  [`KeepTip()`](https://ms609.github.io/TreeTools/reference/DropTip.html).

## Value

A vector specifying an integer, for each tree, which of `tips[-1]` is
most closely related to `tips[1]`.

## See also

Other utility functions:
[`ClusterStrings()`](https://ms609.github.io/TreeSearch/dev/reference/ClusterStrings.md),
[`QACol()`](https://ms609.github.io/TreeSearch/dev/reference/QACol.md),
[`WhenFirstHit()`](https://ms609.github.io/TreeSearch/dev/reference/WhenFirstHit.md)

## Examples

``` r
trees <- inapplicable.trees[["Vinther2008"]]
tips <- c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")
QuartetResolution(trees, tips)
#>  [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 3 2 3 3 2 1 3 3 3 2 3 2 2 2
#> [39] 3 3 2 2 1 2 3 2 1 1 2 1 3 2 3 3 2 1 1 1 3 1 2 1 2 1 3 3 2 1 2 1 2 2 3
```
