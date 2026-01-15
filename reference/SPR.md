# Subtree pruning and rearrangement (SPR)

Perform one SPR rearrangement on a tree

## Usage

``` r
SPR(tree, edgeToBreak = NULL, mergeEdge = NULL)

SPRMoves(tree, edgeToBreak = integer(0))

# S3 method for class 'phylo'
SPRMoves(tree, edgeToBreak = integer(0))

# S3 method for class 'matrix'
SPRMoves(tree, edgeToBreak = integer(0))

SPRSwap(
  parent,
  child,
  nEdge = length(parent),
  nNode = nEdge/2L,
  edgeToBreak = NULL,
  mergeEdge = NULL
)

RootedSPR(tree, edgeToBreak = NULL, mergeEdge = NULL)

RootedSPRSwap(
  parent,
  child,
  nEdge = length(parent),
  nNode = nEdge/2L,
  edgeToBreak = NULL,
  mergeEdge = NULL
)
```

## Arguments

- tree:

  A bifurcating tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), with all nodes
  resolved;

- edgeToBreak:

  the index of an edge to bisect, generated randomly if not specified.

- mergeEdge:

  the index of an edge on which to merge the broken edge.

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

- nNode:

  (optional) Number of nodes.

## Value

This function returns a tree in `phyDat` format that has undergone one
SPR iteration.

[`TBRMoves()`](https://ms609.github.io/TreeSearch/reference/TBR.md)
returns a list of all trees one SPR move away from `tree`, with edges
and nodes in preorder, rooted on the first-labelled tip.

a list containing two elements, corresponding in turn to the rearranged
parent and child parameters

a list containing two elements, corresponding in turn to the rearranged
parent and child parameters

## Details

Equivalent to `kSPR()` in the phangorn package, but faster. Note that
rearrangements that only change the position of the root WILL be
returned by `SPR`. If the position of the root is irrelevant (as in
Fitch parsimony, for example) then this function will occasionally
return a functionally equivalent topology. `RootIrrelevantSPR` will
search tree space more efficiently in these cases. Branch lengths are
not (yet) supported.

All nodes in a tree must be bifurcating;
[ape::collapse.singles](https://rdrr.io/pkg/ape/man/collapse.singles.html)
and [ape::multi2di](https://rdrr.io/pkg/ape/man/multi2di.html) may help.

## Functions

- `SPRSwap()`: faster version that takes and returns parent and child
  parameters

- `RootedSPR()`: Perform SPR rearrangement, retaining position of root

- `RootedSPRSwap()`: faster version that takes and returns parent and
  child parameters

## References

The SPR algorithm is summarized in Felsenstein J (2004). *Inferring
phylogenies*. Sinauer Associates, Sunderland, Massachusetts.

## See also

- `RootedSPR()`: useful when the position of the root node should be
  retained.

Other tree rearrangement functions:
[`NNI()`](https://ms609.github.io/TreeSearch/reference/NNI.md),
[`TBR()`](https://ms609.github.io/TreeSearch/reference/TBR.md)

## Author

Martin R. Smith

## Examples

``` r
{
tree <- ape::rtree(20, br=FALSE)
SPR(tree)
}
#> 
#> Phylogenetic tree with 20 tips and 19 internal nodes.
#> 
#> Tip labels:
#>   t17, t9, t2, t12, t5, t6, ...
#> 
#> Rooted; includes branch length(s).
```
