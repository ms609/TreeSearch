# Minimum and Maximum lengths possible for a character

The smallest and largest length that a phylogenetic character can attain
on any tree.

## Usage

``` r
MinimumLength(x, compress = FALSE)

# S3 method for class 'phyDat'
MinimumLength(x, compress = FALSE)

# S3 method for class 'numeric'
MinimumLength(x, compress = NA)

# S3 method for class 'character'
MinimumLength(x, compress = TRUE)

# S3 method for class 'character'
MaximumLength(x, compress = TRUE)

MinimumSteps(x)

MaximumLength(x, compress = TRUE)

# S3 method for class 'numeric'
MaximumLength(x, compress = NA)
```

## Arguments

- x:

  An object of class `phyDat`; or a string to be coerced to a `phyDat`
  object via
  [`TreeTools::StringToPhyDat()`](https://ms609.github.io/TreeTools/reference/PhyToString.html);
  or an integer vector listing the tokens that may be present at each
  tip along a single character, with each token represented as a binary
  digit; e.g. a value of 11 ( = 2^0 + 2^1 + 2^3) means that the tip may
  have tokens 0, 1 or 3.

  Inapplicable tokens should be denoted with the integer `0` (not 2^0).

- compress:

  Logical specifying whether to retain the compression of a `phyDat`
  object or to return a vector specifying to each individual character,
  decompressed using the dataset's `index` attribute.

## Value

`MinimumLength()` returns a vector of integers specifying the minimum
number of steps that each character must contain.

`MaximumLength()` returns a vector of integers specifying the maximum
number of steps that each character can attain in a parsimonious
reconstruction on a tree. Inapplicable tokens are not yet supported.

## Details

Ambiguous inapplicable states (e.g. `{0, -}`) are currently replaced
with the plain inapplicable token `-`, reflecting the current behaviour
of Morphy.

## See also

Other tree scoring:
[`CharacterLength()`](https://ms609.github.io/TreeSearch/dev/reference/CharacterLength.md),
[`ExpectedLength()`](https://ms609.github.io/TreeSearch/dev/reference/ExpectedLength.md),
[`IWScore()`](https://ms609.github.io/TreeSearch/dev/reference/TreeLength.md),
[`LengthAdded()`](https://ms609.github.io/TreeSearch/dev/reference/LengthAdded.md),
[`MorphyTreeLength()`](https://ms609.github.io/TreeSearch/dev/reference/MorphyTreeLength.md),
[`TaxonInfluence()`](https://ms609.github.io/TreeSearch/dev/reference/TaxonInfluence.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("inapplicable.datasets")
myPhyDat <- inapplicable.phyData[[4]]

# load your own data with
# my.PhyDat <- as.phyDat(read.nexus.data("filepath"))
# or Windows users can select a file interactively using:
# my.PhyDat <- as.phyDat(read.nexus.data(choose.files()))

class(myPhyDat) # phyDat object
#> [1] "phyDat"

# Minimum length of each character in turn
MinimumLength(myPhyDat)
#>   [1] 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 1 1 1 1 1 2 1 1 2 1 1 2 1 1 4 3 1 1 1 1 2 1
#>  [75] 1 1 1 1 1 1 2 4 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [112] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# Collapse duplicate characters, per phyDat compression
MinimumLength(myPhyDat, compress = TRUE)
#>   [1] 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 1 1 1 1 1 2 1 1 2 1 1 2 1 1 4 3 1 1 1 1 2 1
#>  [75] 1 1 1 1 1 1 2 4 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [112] 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# Calculate length of a single character from its textual representation
MinimumLength("-{-1}{-2}{-3}2233")
#> [1] 2
MaximumLength("----0011")
#> [1] 3
```
