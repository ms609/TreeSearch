# Double NNI

Returns the edge parameter of the two trees consistent with the
speficied NNI rearrangement

## Usage

``` r
DoubleNNI(parent, child, edgeToBreak)
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

- edgeToBreak:

  In
  (`Rooted`/`Double`)[`NNI()`](https://ms609.github.io/TreeSearch/reference/NNI.md),
  an optional integer specifying the index of an edge to rearrange,
  generated randomly if not specified. If `-1`, a complete list of all
  trees one step from the input tree will be returned. In
  [`cNNI()`](https://ms609.github.io/TreeSearch/reference/NNI.md), an
  integer from zero to `nEdge(tree) - nTip(tree) - 2`, specifying which
  internal edge to break.

## Value

the `tree[["edge"]]` parameter of the two trees consistent with the
specified rearrangement

## Author

Martin R. Smith
