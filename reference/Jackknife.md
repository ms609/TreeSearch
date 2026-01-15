# Jackknife resampling

Resample trees using Jackknife resampling, i.e. removing a subset of
characters.

## Usage

``` r
Jackknife(
  tree,
  dataset,
  resampleFreq = 2/3,
  InitializeData = PhyDat2Morphy,
  CleanUpData = UnloadMorphy,
  TreeScorer = MorphyLength,
  EdgeSwapper = TBRSwap,
  jackIter = 5000L,
  searchIter = 4000L,
  searchHits = 42L,
  verbosity = 1L,
  ...
)
```

## Arguments

- tree:

  A fully-resolved starting tree in
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html) format, with the
  desired outgroup. Edge lengths are not supported and will be removed.

- dataset:

  a dataset in the format required by `TreeScorer()`.

- resampleFreq:

  Double between 0 and 1 stating proportion of characters to resample.

- InitializeData:

  Function that sets up data object to prepare for tree search. The
  function will be passed the `dataset` parameter. Its return value will
  be passed to `TreeScorer()` and `CleanUpData()`.

- CleanUpData:

  Function to destroy data object on function exit. The function will be
  passed the value returned by `InitializeData()`.

- TreeScorer:

  function to score a given tree. The function will be passed three
  parameters, corresponding to the `parent` and `child` entries of a
  tree's edge list, and a dataset.

- EdgeSwapper:

  a function that rearranges a parent and child vector, and returns a
  list with modified vectors; for example
  [`SPRSwap()`](https://ms609.github.io/TreeSearch/reference/SPR.md).

- jackIter:

  Integer specifying number of jackknife iterations to conduct.

- searchIter:

  Integer specifying maximum rearrangements to perform on each bootstrap
  or ratchet iteration. To override this value for a single swapper
  function, set e.g. `attr(SwapperFunction, "searchIter") <- 99`

- searchHits:

  Integer specifying maximum times to hit best score before terminating
  a tree search within a ratchet iteration. To override this value for a
  single swapper function, set e.g.
  `attr(SwapperFunction, "searchHits") <- 99`

- verbosity:

  Numeric specifying level of detail to display in console: larger
  numbers provide more verbose feedback to the user.

- ...:

  further parameters to send to `TreeScorer()`

## Value

`Jackknife()` returns a list of trees recovered after jackknife
iterations.

## Details

The function assumes that `InitializeData()` will return a morphy
object; if this doesn't hold for you, post a [GitHub
issue](https://github.com/ms609/TreeSearch/issues/new/) or e-mail the
maintainer.

## See also

- [`Resample()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md):
  Jackknife resampling for non-custom searches performed using
  [`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md).

- [`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md):
  Label nodes of a tree with jackknife supports.

Other split support functions:
[`ConcordanceTable()`](https://ms609.github.io/TreeSearch/reference/ConcordanceTable.md),
[`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md),
[`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md),
[`MostContradictedFreq()`](https://ms609.github.io/TreeSearch/reference/MostContradictedFreq.md),
[`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md),
[`SiteConcordance`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)

Other custom search functions:
[`EdgeListSearch()`](https://ms609.github.io/TreeSearch/reference/TreeSearch.md),
[`MorphyBootstrap()`](https://ms609.github.io/TreeSearch/reference/Ratchet.md),
[`SuccessiveApproximations()`](https://ms609.github.io/TreeSearch/reference/SuccessiveApproximations.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)
