# Non-duplicate root

Identify, for each edge, whether it denotes a different partition from
the root edge. The first edge of the input tree must be a root edge;
this can be accomplished using
[`Preorder()`](https://ms609.github.io/TreeTools/reference/Reorder.html).

## Usage

``` r
.NonDuplicateRoot(parent, child, nEdge = length(parent))
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

  (optional) integer specifying the number of edges of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `dim(tree$edge)[1]`

## Value

`.NonDuplicateRoot()` returns a logical vector of length `nEdge`,
specifying `TRUE` unless an edge identifies the same partition as the
root edge.

## Details

This function is a copy of a deprecated ancestor in TreeTools; see
[\#32](https://github.com/ms609/TreeTools/issues/32).

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- TreeTools::Preorder(TreeTools::BalancedTree(8))
edge <- tree$edge
parent <- edge[, 1]
child <- edge[, 2]

which(!.NonDuplicateRoot(parent, child))
#> [1] 1
```
