# Empirically counted profiles for small trees

The base 2 logarithm of the number of trees containing *s* steps,
calculated by scoring a character on each *n*-leaf tree.

## Usage

``` r
profiles
```

## Format

A list with the structure
`profiles[[number of leaves]][[number of tokens]][[tokens in smallest split]]`
The list entry returns a named numeric vector; each entry lists
log2(proportion of *n*-leaf trees with *s* or fewer steps for this
character).

## See also

Other profile parsimony functions:
[`Carter1()`](https://ms609.github.io/TreeSearch/dev/reference/Carter1.md),
[`PrepareDataProfile()`](https://ms609.github.io/TreeSearch/dev/reference/PrepareDataProfile.md),
[`StepInformation()`](https://ms609.github.io/TreeSearch/dev/reference/StepInformation.md),
[`WithOneExtraStep()`](https://ms609.github.io/TreeSearch/dev/reference/WithOneExtraStep.md)

## Examples

``` r
data(profiles)

# Load profile for a character of the structure 0 0 0 1 1 1 1 1
profile3.5 <- profiles[[8]][[2]][[3]]

# Number of trees with _s_ or fewer steps on that character
TreeTools::NUnrooted(8) * 2 ^ profile3.5
#>     1     2     3     4 
#>   225  2475  8019 10395 
```
