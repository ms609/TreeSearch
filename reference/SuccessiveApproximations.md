# Tree search using successive approximations

Searches for a tree that is optimal under the Successive Approximations
criterion (Farris 1969) .

## Usage

``` r
SuccessiveApproximations(
  tree,
  dataset,
  outgroup = NULL,
  k = 3,
  maxSuccIter = 20,
  ratchetHits = 100,
  searchHits = 50,
  searchIter = 500,
  ratchetIter = 5000,
  verbosity = 0,
  suboptimal = 0.1
)

SuccessiveWeights(tree, dataset)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree. Perhaps load into R
  using
  [`ReadAsPhyDat()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html).
  Additive (ordered) characters can be handled using
  [`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html).

- outgroup:

  if not NULL, taxa on which the tree should be rooted

- k:

  Constant for successive approximations, see Farris 1969 p. 379

- maxSuccIter:

  maximum iterations of successive approximation

- ratchetHits:

  maximum hits for parsimony ratchet

- searchHits:

  maximum hits in tree search

- searchIter:

  maximum iterations in tree search

- ratchetIter:

  maximum iterations of parsimony ratchet

- verbosity:

  Integer specifying level of messaging; higher values give more
  detailed commentary on search progress. Set to `0` to run silently.

- suboptimal:

  retain trees that are this proportion less optimal than the optimal
  tree

## Value

`SuccessiveApproximations()` returns a list of class `multiPhylo`
containing optimal (and slightly suboptimal, if suboptimal \> 0) trees.

`SuccessiveWeights()` returns the score of a tree, given the weighting
instructions specified in the attributes of the dataset.

## References

Farris JS (1969). “A successive approximations approach to character
weighting.” *Systematic Biology*, **18**(4), 374–385.
[doi:10.2307/2412182](https://doi.org/10.2307/2412182) .

## See also

Other custom search functions:
[`EdgeListSearch()`](https://ms609.github.io/TreeSearch/reference/TreeSearch.md),
[`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md),
[`MorphyBootstrap()`](https://ms609.github.io/TreeSearch/reference/Ratchet.md)
