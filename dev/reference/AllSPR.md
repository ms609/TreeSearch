# All SPR trees

All SPR trees

## Usage

``` r
AllSPR(parent, child, nEdge, notDuplicateRoot, edgeToBreak)
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

- nEdge:

  integer specifying the number of edges of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `dim(tree$edge)[1]`

- notDuplicateRoot:

  logical vector of length `nEdge`, specifying for each edge whether it
  is the second edge leading to the root (in which case its breaking
  will be equivalent to breaking the other root edge... except insofar
  as it moves the position of the root.)

- edgeToBreak:

  (optional) integer specifying the index of an edge to bisect/prune,
  generated randomly if not specified. Alternatively, set to `-1` to
  return a complete list of all trees one step from the input tree.

## Value

`AllSPR()` returns a list of edge matrices for all trees one SPR
rearrangement from the starting tree

## Author

Martin R. Smith
