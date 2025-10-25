# Parsimony score of random postorder tree

Parsimony score of random postorder tree

## Usage

``` r
RandomTreeScore(morphyObj)
```

## Arguments

- morphyObj:

  Object of class `morphy`, perhaps created with
  [`PhyDat2Morphy()`](https://ms609.github.io/TreeSearch/reference/PhyDat2Morphy.md).

## Value

`RandomTreeScore()` returns the parsimony score of a random tree for the
given Morphy object.

## Examples

``` r
tokens <- matrix(c(
  0, "-", "-", 1, 1, 2,
  0, 1, 0, 1, 2, 2,
  0, "-", "-", 0, 0, 0), byrow = TRUE, nrow = 3L,
  dimnames = list(letters[1:3], NULL))
pd <- TreeTools::MatrixToPhyDat(tokens)
morphyObj <- PhyDat2Morphy(pd)

RandomTreeScore(morphyObj)
#> [1] 4

morphyObj <- UnloadMorphy(morphyObj)
```
