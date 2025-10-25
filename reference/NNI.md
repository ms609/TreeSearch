# Nearest neighbour interchange (NNI)

`NNI()`performs a single iteration of the nearest-neighbour interchange
algorithm; `RootedNNI()` retains the position of the root. These
functions are based on equivalents in the phangorn package. `cNNI()` is
an equivalent function coded in C, that runs much faster.

## Usage

``` r
NNI(tree, edgeToBreak = NULL)

cNNI(tree, edgeToBreak = NULL, whichSwitch = NULL)

NNISwap(parent, child, nTips = (length(parent)/2L) + 1L, edgeToBreak = NULL)

RootedNNI(tree, edgeToBreak = NULL)

RootedNNISwap(
  parent,
  child,
  nTips = (length(parent)/2L) + 1L,
  edgeToBreak = NULL
)
```

## Arguments

- tree:

  A tree of class `phylo`. For `cNNI()`, this must be a binary tree
  rooted on a single leaf, whose root node is the lowest numbered
  internal node.

- edgeToBreak:

  In (`Rooted`/`Double`)`NNI()`, an optional integer specifying the
  index of an edge to rearrange, generated randomly if not specified. If
  `-1`, a complete list of all trees one step from the input tree will
  be returned. In `cNNI()`, an integer from zero to
  `nEdge(tree) - nTip(tree) - 2`, specifying which internal edge to
  break.

- whichSwitch:

  Integer from zero to one, specifying which way to re-build the broken
  internal edge.

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

- child:

  Integer vector corresponding to the second column of the edge matrix
  of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `tree[["edge"]][, 2]`.

- nTips:

  (optional) Number of tips.

## Value

Returns a tree with class `phylo` (if `returnAll = FALSE`) or a set of
trees, with class `multiPhylo` (if `returnAll = TRUE`).

`cNNI()` returns a tree of class `phylo`, rooted on the same leaf, on
which the specified rearrangement has been conducted.

`NNISwap()` returns a list containing two elements, corresponding in
turn to the rearranged parent and child parameters.

a list containing two elements, corresponding in turn to the rearranged
parent and child parameters

## Details

Branch lengths are not supported.

All nodes in a tree must be bifurcating;
[`ape::collapse.singles()`](https://rdrr.io/pkg/ape/man/collapse.singles.html)
and [`ape::multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html) may
help.

## Functions

- `NNISwap()`: faster version that takes and returns parent and child
  parameters

- `RootedNNI()`: Perform NNI rearrangement, retaining position of root

- `RootedNNISwap()`: faster version that takes and returns parent and
  child parameters

## References

The algorithm is summarized in Felsenstein J (2004). *Inferring
phylogenies*. Sinauer Associates, Sunderland, Massachusetts.

## See also

Other tree rearrangement functions:
[`SPR()`](https://ms609.github.io/TreeSearch/reference/SPR.md),
[`TBR()`](https://ms609.github.io/TreeSearch/reference/TBR.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- TreeTools::BalancedTree(8)
# A random rearrangement
NNI(tree)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
cNNI(tree)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.

# All trees one NNI rearrangement away
NNI(tree, edgeToBreak = -1)
#> 12 phylogenetic trees

# Manual random sampling
cNNI(tree, sample.int(14 - 8 - 1, 1), sample.int(2, 1))
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.

# A specified rearrangement
cNNI(tree, 0, 0)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.

# If a tree may not be binary, collapse nodes with
tree <- TreeTools::MakeTreeBinary(tree)

# If a tree may be improperly rooted, use
tree <- TreeTools::RootTree(tree, 1)

# If a tree may exhibit unusual node ordering, this can be addressed with
tree <- TreeTools::Preorder(tree)
```
