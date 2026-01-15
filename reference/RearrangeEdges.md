# Rearrange edges of a phylogenetic tree

`RearrangeEdges()` performs the specified edge rearrangement on a matrix
that corresponds to the edges of a phylogenetic tree, returning the
score of the new tree. Will generally be called from within a tree
search function.

## Usage

``` r
RearrangeEdges(
  parent,
  child,
  dataset,
  TreeScorer = MorphyLength,
  EdgeSwapper,
  scoreToBeat = TreeScorer(parent, child, dataset, ...),
  iter = "?",
  hits = 0L,
  verbosity = 0L,
  ...
)
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

- dataset:

  Third argument to pass to `TreeScorer`.

- TreeScorer:

  function to score a given tree. The function will be passed three
  parameters, corresponding to the `parent` and `child` entries of a
  tree's edge list, and a dataset.

- EdgeSwapper:

  a function that rearranges a parent and child vector, and returns a
  list with modified vectors; for example
  [`SPRSwap()`](https://ms609.github.io/TreeSearch/reference/SPR.md).

- scoreToBeat:

  Double giving score of input tree.

- iter:

  iteration number of calling function, for reporting to user only.

- hits:

  Integer giving number of times the input tree has already been hit.

- verbosity:

  Numeric specifying level of detail to display in console: larger
  numbers provide more verbose feedback to the user.

- ...:

  further arguments to pass to `TreeScorer()`, e.g. `dataset = `.

## Value

This function returns a list with two to four elements, corresponding to
a binary tree: - 1. Integer vector listing the parent node of each
edge; - 2. Integer vector listing the child node of each edge; - 3.
Score of the tree; - 4. Number of times that score has been hit.

## Details

`RearrangeTree()` performs one tree rearrangement of a specified type,
and returns the score of the tree (with the given dataset). It also
reports the number of times that this score was hit in the current
function call.

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("Lobo", package="TreeTools")
tree <- TreeTools::NJTree(Lobo.phy)
edge <- tree$edge
parent <- edge[, 1]
child <- edge[, 2]
dataset <- PhyDat2Morphy(Lobo.phy)
RearrangeEdges(parent, child, dataset, EdgeSwapper = RootedNNISwap)
#> [[1]]
#>  [1] 49 49 50 50 51 52 53 54 55 56 57 57 56 58 59 59 58 60 61 61 60 62 63 63 64
#> [26] 64 65 65 66 67 67 66 68 68 62 69 70 71 72 72 71 70 73 73 69 55 74 74 75 76
#> [51] 76 75 54 53 77 77 78 79 79 78 80 81 82 82 81 80 83 83 84 85 86 86 85 87 87
#> [76] 84 52 88 88 51 89 90 91 92 92 91 93 93 94 94 95 95 90 89
#> 
#> [[2]]
#>  [1]  1 50  2 51 52 53 54 55 56 57  3  4 58 59 31 32 60 61 34 35 62 63 36 64 43
#> [26] 65 44 66 67 45 46 68 47 48 69 70 71 72 37 41 42 73 38 39 40 74  5 75 76 28
#> [51] 29 30  6 77  9 78 79 10 11 80 81 82 12 13 22 83 14 84 85 86 15 18 87 16 17
#> [76] 19 88  7  8 89 90 91 92 20 21 93 24 94 25 95 26 27 23 33
#> 
#> [[3]]
#> [1] 232
#> 
#> [[4]]
#> [1] 0
#> 
# Remember to free memory:
dataset <- UnloadMorphy(dataset)
```
