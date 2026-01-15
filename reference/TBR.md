# Tree bisection and reconnection (TBR)

`TBR` performs a single random TBR iteration.

## Usage

``` r
TBR(tree, edgeToBreak = NULL, mergeEdges = NULL)

TBRMoves(tree, edgeToBreak = integer(0))

# S3 method for class 'phylo'
TBRMoves(tree, edgeToBreak = integer(0))

# S3 method for class 'matrix'
TBRMoves(tree, edgeToBreak = integer(0))

TBRSwap(
  parent,
  child,
  nEdge = length(parent),
  edgeToBreak = NULL,
  mergeEdges = NULL
)

RootedTBR(tree, edgeToBreak = NULL, mergeEdges = NULL)

RootedTBRSwap(
  parent,
  child,
  nEdge = length(parent),
  edgeToBreak = NULL,
  mergeEdges = NULL
)
```

## Arguments

- tree:

  A bifurcating tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), with all nodes
  resolved;

- edgeToBreak:

  (optional) integer specifying the index of an edge to bisect/prune,
  generated randomly if not specified. Alternatively, set to `-1` to
  return a complete list of all trees one step from the input tree.

- mergeEdges:

  (optional) vector of length 1 or 2, listing edge(s) to be joined: In
  SPR, this is where the pruned subtree will be reconnected. In TBR,
  these edges will be reconnected (so must be on opposite sides of
  `edgeToBreak`); if only a single edge is specified, the second will be
  chosen at random

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

  (optional) Number of edges.

## Value

`TBR()` returns a tree in `phyDat` format that has undergone one TBR
iteration.

`TBRMoves()` returns a `multiPhylo` object listing all trees one TBR
move away from `tree`, with edges and nodes in preorder, rooted on the
first-labelled tip.

`TBRSwap()` returns a list containing two elements corresponding to the
rearranged `parent` and `child` parameters.

## Details

Branch lengths are not (yet) supported.

All nodes in a tree must be bifurcating;
[ape::collapse.singles](https://rdrr.io/pkg/ape/man/collapse.singles.html)
and [ape::multi2di](https://rdrr.io/pkg/ape/man/multi2di.html) may help.

## Functions

- `TBRSwap()`: faster version that takes and returns parent and child
  parameters

- `RootedTBR()`: Perform TBR rearrangement, retaining position of root

- `RootedTBRSwap()`: faster version that takes and returns parent and
  child parameters

## References

The TBR algorithm is summarized in Felsenstein J (2004). *Inferring
phylogenies*. Sinauer Associates, Sunderland, Massachusetts.

## See also

`RootedTBR()`: useful when the position of the root node should be
retained.

Other tree rearrangement functions:
[`NNI()`](https://ms609.github.io/TreeSearch/reference/NNI.md),
[`SPR()`](https://ms609.github.io/TreeSearch/reference/SPR.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
library("ape")
tree <- rtree(20, br=NULL)
TBR(tree)
#> 
#> Phylogenetic tree with 20 tips and 19 internal nodes.
#> 
#> Tip labels:
#>   t7, t5, t20, t15, t3, t17, ...
#> 
#> Rooted; no branch length.
```
