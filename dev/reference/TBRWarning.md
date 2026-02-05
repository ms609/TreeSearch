# TBR Warning Print a warning and return given tree

TBR Warning Print a warning and return given tree

## Usage

``` r
SPRWarning(parent, child, error)

TBRWarning(parent, child, error)
```

## Arguments

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

- child:

  Integer vector corresponding to the second column of the edge matrix
  of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `tree[["edge"]][, 2]`.

- error:

  error message to report

## Value

A list with the entries `parent`, `child`.

## Functions

- `SPRWarning()`: for SPR rearrangements

## Author

Martin R. Smith

## Examples

``` r
suppressWarnings(TBRWarning(0, 0, "Message text")) # will trigger warning
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 0
#> 
```
