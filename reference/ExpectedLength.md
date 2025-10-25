# Expected length

For a given dataset and tree topology, `ExpectedLength()` estimates the
length expected if the states of each character are shuffled randomly
across the leaves.

## Usage

``` r
ExpectedLength(dataset, tree, nRelabel = 1000, compress = FALSE)
```

## Arguments

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree. Perhaps load into R
  using
  [`ReadAsPhyDat()`](https://ms609.github.io/TreeTools/reference/ReadCharacters.html).
  Additive (ordered) characters can be handled using
  [`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.html).

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- nRelabel:

  Integer specifying how many times to relabel leaves when estimating
  null tree length for RHI calculation. Steell et al. (2025) recommend
  1000, but suggest that 100 may suffice. If zero (the default), the RHI
  is not calculated.

- compress:

  Logical specifying whether to retain the compression of a `phyDat`
  object or to return a vector specifying to each individual character,
  decompressed using the dataset's `index` attribute.

## Value

`ExpectedLength()` returns a numeric vector stating the median length of
each character in `dataset` on `tree` after `nRelabel` random
relabelling of leaves.

## References

Steell EM, Hsiang AY, Field DJ (2025). “Revealing Patterns of Homoplasy
in Discrete Phylogenetic Datasets with a Cross-Comparable Index.”
*Zoological Journal of the Linnean Society*, **204**(1), zlaf024.
[doi:10.1093/zoolinnean/zlaf024](https://doi.org/10.1093/zoolinnean/zlaf024)
.

## See also

Other tree scoring:
[`CharacterLength()`](https://ms609.github.io/TreeSearch/reference/CharacterLength.md),
[`IWScore()`](https://ms609.github.io/TreeSearch/reference/TreeLength.md),
[`LengthAdded()`](https://ms609.github.io/TreeSearch/reference/LengthAdded.md),
[`MinimumLength()`](https://ms609.github.io/TreeSearch/reference/MinimumLength.md),
[`MorphyTreeLength()`](https://ms609.github.io/TreeSearch/reference/MorphyTreeLength.md),
[`TaxonInfluence()`](https://ms609.github.io/TreeSearch/reference/TaxonInfluence.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)
