# Combine two edge matrices

Combine two edge matrices

## Usage

``` r
.CombineResults(x, y, stage)

.ReplaceResults(old, new, stage)
```

## Arguments

- x, y:

  3D arrays, each slice containing an edge matrix from a tree of class
  `phylo`. `x` should not contain duplicates.

- stage:

  Integer specifying element of `firstHit` in which new hits should be
  recorded.

- old:

  old array of edge matrices with `firstHit` attribute.

- new:

  new array of edge matrices.

## Value

A single 3D array containing each unique edge matrix from (`x` and) `y`,
with a `firstHit` attribute as documented in
[`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md).

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)
