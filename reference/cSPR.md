# `cSPR()` expects a tree rooted on a single tip.

`cSPR()` expects a tree rooted on a single tip.

## Usage

``` r
cSPR(tree, whichMove = NULL)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- whichMove:

  Integer specifying which SPR move index to perform.

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tree <- TreeTools::BalancedTree(8)

# Tree must be rooted on leaf
tree <- TreeTools::RootTree(tree, 1)

# Random rearrangement
cSPR(tree)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.

# Specific rearrangement
cSPR(tree, 9)
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#> 
#> Rooted; no branch length.
```
