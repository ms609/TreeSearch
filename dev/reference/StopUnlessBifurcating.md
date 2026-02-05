# Check that all nodes in a tree are bifurcating.

Check that all nodes in a tree are bifurcating.

## Usage

``` r
StopUnlessBifurcating(parent)
```

## Arguments

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

## Value

Returns `NULL`, but will `stop` with an error message if a tree does not
appear to be bifurcating.

## Author

Martin R. Smith
