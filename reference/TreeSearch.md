# Search for most parsimonious trees

Run standard search algorithms (NNI, SPR or TBR) to search for a more
parsimonious tree.

For detailed documentation of the "TreeSearch" package, including full
instructions for loading phylogenetic data into R and initiating and
configuring tree search, see the [package
documentation](https://ms609.github.io/TreeSearch/).

## Usage

``` r
EdgeListSearch(
  edgeList,
  dataset,
  TreeScorer = MorphyLength,
  EdgeSwapper = RootedTBRSwap,
  maxIter = 100,
  maxHits = 20,
  bestScore = NULL,
  stopAtScore = NULL,
  stopAtPeak = FALSE,
  stopAtPlateau = 0L,
  verbosity = 1L,
  ...
)

TreeSearch(
  tree,
  dataset,
  InitializeData = PhyDat2Morphy,
  CleanUpData = UnloadMorphy,
  TreeScorer = MorphyLength,
  EdgeSwapper = RootedTBRSwap,
  maxIter = 100L,
  maxHits = 20L,
  stopAtPeak = FALSE,
  stopAtPlateau = 0L,
  verbosity = 1L,
  ...
)

IWTreeSearch(...)

EmptyPhyDat(tree)

DoNothing(x = NULL)
```

## Arguments

- edgeList:

  a list containing the following:

  - vector of integers corresponding to the parent of each edge in turn

  - vector of integers corresponding to the child of each edge in turn

  - (optionally) score of the tree

  - (optionally, if score provided) number of times this score has been
    hit

- dataset:

  Data in format required by `InitializeData`.

- TreeScorer:

  function to score a given tree. The function will be passed three
  parameters, corresponding to the `parent` and `child` entries of a
  tree's edge list, and a dataset.

- EdgeSwapper:

  a function that rearranges a parent and child vector, and returns a
  list with modified vectors; for example
  [`SPRSwap()`](https://ms609.github.io/TreeSearch/reference/SPR.md).

- maxIter:

  Numeric specifying maximum number of iterations to perform before
  abandoning the search.

- maxHits:

  Numeric specifying maximum times to hit the best pscore before
  abandoning the search.

- stopAtPeak:

  Logical specifying whether to terminate search once a subsequent
  iteration recovers a sub-optimal score. Will be overridden if a passed
  function has an attribute `stopAtPeak` set by
  `attr(FunctionName, "stopAtPeak") <- TRUE`.

- stopAtPlateau:

  Integer. If \> 0, tree search will terminate if the score has not
  improved after `stopAtPlateau` iterations. Will be overridden if a
  passed function has an attribute `stopAtPlateau` set by
  `attr(FunctionName, "stopAtPlateau") <- TRUE`.

- verbosity:

  Numeric specifying level of detail to display in console: larger
  numbers provide more verbose feedback to the user.

- ...:

  further arguments to pass to `TreeScorer()`, e.g. `dataset = `.

- tree:

  A fully-resolved starting tree in
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html) format, with the
  desired outgroup. Edge lengths are not supported and will be removed.

- InitializeData:

  Function that sets up data object to prepare for tree search. The
  function will be passed the `dataset` parameter. Its return value will
  be passed to `TreeScorer()` and `CleanUpData()`.

- CleanUpData:

  Function to destroy data object on function exit. The function will be
  passed the value returned by `InitializeData()`.

## Value

`TreeSearch()` returns a tree, with an attribute `pscore` conveying its
parsimony score. \#" Note that the parsimony score will be inherited
from the tree"s attributes, which is only valid if it was generated
using the same `data` that is passed here.

`EmptyPhyDat()` returns a `phyDat` object comprising a single null
character, coded with state zero for every leaf in `tree`.

## Functions

- `EdgeListSearch()`: Tree search from edge lists

## See also

- [`Fitch`](https://ms609.github.io/TreeSearch/reference/TreeLength.md),
  calculates parsimony score;

- [`RootedNNI`](https://ms609.github.io/TreeSearch/reference/NNI.md),
  conducts tree rearrangements;

- [`Ratchet`](https://ms609.github.io/TreeSearch/reference/Ratchet.md),
  alternative heuristic, useful to escape local optima.

Other custom search functions:
[`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md),
[`MorphyBootstrap()`](https://ms609.github.io/TreeSearch/reference/Ratchet.md),
[`SuccessiveApproximations()`](https://ms609.github.io/TreeSearch/reference/SuccessiveApproximations.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("Lobo", package="TreeTools")
njtree <- TreeTools::NJTree(Lobo.phy)

## Only run examples in interactive R sessions
if (interactive()) {
  TreeSearch(njtree, Lobo.phy, maxIter = 20, EdgeSwapper = NNISwap)
  TreeSearch(njtree, Lobo.phy, maxIter = 20, EdgeSwapper = RootedSPRSwap)
  TreeSearch(njtree, Lobo.phy, maxIter = 20, EdgeSwapper = TBRSwap)
}
```
