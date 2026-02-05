# Character length

Homoplasy length of each character in a dataset on a specified tree.

## Usage

``` r
CharacterLength(tree, dataset, compress = FALSE)

FastCharacterLength(tree, dataset)
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

- compress:

  Logical specifying whether to retain the compression of a `phyDat`
  object or to return a vector specifying to each individual character,
  decompressed using the dataset's `index` attribute.

## Value

`CharacterLength()` returns a vector listing the contribution of each
character to tree score, according to the algorithm of Brazeau et al.
(2019) .

## Functions

- `FastCharacterLength()`: Do not perform checks. Use with care: may
  cause erroneous results or software crash if variables are in the
  incorrect format.

## References

Brazeau MD, Guillerme T, Smith MR (2019). “An algorithm for
morphological phylogenetic analysis with inapplicable data.” *Systematic
Biology*, **68**(4), 619–631.
[doi:10.1093/sysbio/syy083](https://doi.org/10.1093/sysbio/syy083) .

## See also

Other tree scoring:
[`ExpectedLength()`](https://ms609.github.io/TreeSearch/dev/reference/ExpectedLength.md),
[`IWScore()`](https://ms609.github.io/TreeSearch/dev/reference/TreeLength.md),
[`LengthAdded()`](https://ms609.github.io/TreeSearch/dev/reference/LengthAdded.md),
[`MinimumLength()`](https://ms609.github.io/TreeSearch/dev/reference/MinimumLength.md),
[`MorphyTreeLength()`](https://ms609.github.io/TreeSearch/dev/reference/MorphyTreeLength.md),
[`TaxonInfluence()`](https://ms609.github.io/TreeSearch/dev/reference/TaxonInfluence.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("inapplicable.datasets")
dataset <- inapplicable.phyData[[12]]
tree <- TreeTools::NJTree(dataset)
CharacterLength(tree, dataset)
#>   [1]  1  1  2  2  3  4  4  5  1  6  1  1  2  1  2  4  2  5  2  4  1  1  2  1  1
#>  [26]  2  2  2  5  1  1  1  2  1  4  3  6  5  5  1  2  1  2  1  1  1  1  1  1  1
#>  [51]  2  1  1  8  3  7  3  3  2  5  2  1  2  2  2  9 10  5  5  2  4  3  5  1  1
#>  [76]  5  6  3  5  5  3  3  2  7  2  1  9  5  7  5  1  7  1  5  2  2  8  1  2  2
#> [101]  1  1  3  2  3  4 12  6  8  4  5  8  7  1  1  2  2  1  1  4  6  2  3  6  6
#> [126]  1  1  1  1  2  1  1  1  2  1  1  6
CharacterLength(tree, dataset, compress = TRUE)
#>   [1]  1  1  2  2  3  4  4  5  1  6  1  1  2  1  2  4  2  5  2  4  1  2  1  1  2
#>  [26]  2  2  5  1  1  1  2  4  3  6  5  5  1  2  2  2  1  8  3  7  3  3  2  5  2
#>  [51]  1  2  2  9 10  5  5  2  4  3  5  1  1  5  6  3  5  5  3  3  2  7  2  1  9
#>  [76]  5  7  5  1  7  1  5  2  2  8  1  2  2  1  3  2  3  4 12  6  8  4  5  8  7
#> [101]  1  2  2  1  4  6  2  3  6  6  1  1  1  1  2  1  2  6
```
