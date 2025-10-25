# Information content of a character known to contain *e* steps

`StepInformation()` calculates the phylogenetic information content of a
character `char` when *e* extra steps are present, for all possible
values of *e*.

## Usage

``` r
StepInformation(char, ambiguousTokens = c("-", "?"))
```

## Arguments

- char:

  Vector of tokens listing states for the character in question.

- ambiguousTokens:

  Vector specifying which tokens, if any, correspond to the ambiguous
  token (`?`).

## Value

`StepInformation()` returns a numeric vector detailing the amount of
phylogenetic information (in bits) associated with the character when 0,
1, 2… extra steps are present. The vector is named with the total number
of steps associated with each entry in the vector: for example, a
character with three observed tokens must exhibit two steps, so the
first entry (zero extra steps) is named `2` (two steps observed).

## Details

Calculates the number of trees consistent with the character having *e*
extra steps, where *e* ranges from its minimum possible value (i.e.
number of different tokens minus one) to its maximum.

## See also

Other profile parsimony functions:
[`Carter1()`](https://ms609.github.io/TreeSearch/reference/Carter1.md),
[`PrepareDataProfile()`](https://ms609.github.io/TreeSearch/reference/PrepareDataProfile.md),
[`WithOneExtraStep()`](https://ms609.github.io/TreeSearch/reference/WithOneExtraStep.md),
[`profiles`](https://ms609.github.io/TreeSearch/reference/profiles.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
character <- rep(c(0:3, "?", "-"), c(8, 5, 1, 1, 2, 2))
StepInformation(character)
#>         3         4         5         6         7 
#> 9.9203529 5.5280354 2.5784492 0.7618403 0.0000000 
```
